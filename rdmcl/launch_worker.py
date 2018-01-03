#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Jul 17 2017 

"""
DESCRIPTION OF PROGRAM
"""

import pandas as pd
import sqlite3
from multiprocessing import Lock, cpu_count
from buddysuite import buddy_resources as br
from buddysuite import AlignBuddy as Alb
from buddysuite import SeqBuddy as Sb
import os
import sys
import re
import time
from random import random
from collections import OrderedDict
import shutil
import traceback

# My packages
try:
    from . import helpers
    from . import rdmcl
except ImportError:
    import helpers
    import rdmcl


# Globals
WORKERLOCK = Lock()
VERSION = helpers.VERSION
VERSION.name = "launch_worker"

# Set global precision levels
pd.set_option("display.precision", 12)


class Worker(object):
    def __init__(self, location, heartrate=60, max_wait=600, dead_thread_wait=120, cpus=cpu_count(),
                 job_size_coff=300, log=False, quiet=False):
        self.working_dir = os.path.abspath(location)
        self.wrkdb_path = os.path.join(self.working_dir, "work_db.sqlite")
        self.hbdb_path = os.path.join(self.working_dir, "heartbeat_db.sqlite")
        self.output = os.path.join(self.working_dir, ".worker_output")
        os.makedirs(self.output, exist_ok=True)

        self.heartrate = heartrate
        self.heartbeat = rdmcl.HeartBeat(self.hbdb_path, self.heartrate, thread_type="worker")
        self.max_wait = max_wait
        self.dead_thread_wait = dead_thread_wait
        self.cpus = cpus - 1 if cpus > 1 else 1
        self.job_size_coff = job_size_coff
        self.worker_file = ""
        self.data_file = ""
        self.start_time = time.time()
        self.split_time = 0
        self.idle = 1
        self.running = 1
        self.last_heartbeat_from_master = 0
        self.subjob_num = 1
        self.num_subjobs = 1
        self.job_id_hash = None
        self.printer = br.DynamicPrint(quiet=quiet)
        if log:
            self.printer._writer = _write(self.printer)

    def start(self):
        self.split_time = time.time()
        self.start_time = time.time()

        self.heartbeat.start()
        self.worker_file = os.path.join(self.working_dir, "Worker_%s" % self.heartbeat.id)
        with open(self.worker_file, "w") as ofile:
            ofile.write("To terminate this Worker, simply delete this file.")

        self.data_file = os.path.join(self.working_dir, ".Worker_%s.dat" % self.heartbeat.id)
        open(self.data_file, "w").close()

        helpers.dummy_func()

        self.last_heartbeat_from_master = time.time()
        self.printer.write("Starting Worker_%s" % self.heartbeat.id)
        self.printer.new_line(1)

        idle_countdown = 1
        while os.path.isfile(self.worker_file):
            idle = round(100 * self.idle / (self.idle + self.running), 2)
            if not idle_countdown:
                self.printer.write("Idle %s%%" % idle)
                idle_countdown = 5

            # Make sure there are some masters still kicking around
            self.check_masters(idle)

            # Check for and clean up dead threads and orphaned jobs every twentieth(ish) time through
            rand_check = random()
            if rand_check > 0.95:
                self.clean_dead_threads()

            # Fetch a job from the queue
            data = self.fetch_queue_job()
            if data:
                full_name, psipred_dir, align_m, align_p, trimal, gap_open, gap_extend = data
                subjob_num, num_subjobs, id_hash = [1, 1, full_name] if len(full_name.split("_")) == 1 \
                    else full_name.split("_")
                subjob_num = int(subjob_num)
                num_subjobs = int(num_subjobs)
                self.printer.write("Running %s" % full_name)
            else:
                time.sleep(random() * self.idle_workers())  # Pause for some time relative to num idle workers
                idle_countdown -= 1
                self.idle += time.time() - self.split_time
                self.split_time = time.time()
                continue

            try:
                idle_countdown = 1
                seqbuddy = Sb.SeqBuddy("%s/%s.seqs" % (self.output, id_hash), in_format="fasta")

                # Prepare alignment
                if len(seqbuddy) == 1:
                    raise ValueError("Queued job of size 1 encountered: %s" % id_hash)
                else:
                    if num_subjobs == 1:
                        self.printer.write("Creating MSA (%s seqs)" % len(seqbuddy))
                        alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), align_m,
                                                     params=align_p, quiet=True)
                    else:
                        self.printer.write("Reading MSA (%s seqs)" % len(seqbuddy))
                        alignment = Alb.AlignBuddy(os.path.join(self.output, "%s.aln" % id_hash))

                # Prepare psipred dataframes
                psipred_dfs = self.prepare_psipred_dfs(seqbuddy, psipred_dir)

                if num_subjobs == 1:  # This is starting a full job from scratch, not a sub-job
                    # Need to specify what columns the PsiPred files map to now that there are gaps.
                    psipred_dfs = rdmcl.update_psipred(alignment, psipred_dfs, "msa")

                    # TrimAl
                    self.printer.write("Trimal (%s seqs)" % len(seqbuddy))
                    alignment = rdmcl.trimal(seqbuddy, trimal, alignment)

                    with helpers.ExclusiveConnect(os.path.join(self.output, "write.lock"), max_lock=0):
                        # Place these write commands in ExclusiveConnect to ensure a writing lock
                        if not os.path.isfile(os.path.join(self.output, "%s.aln" % id_hash)):
                            alignment.write(os.path.join(self.output, "%s.aln" % id_hash), out_format="fasta")

                    # Re-update PsiPred files now that some columns, possibly including non-gap characters, are removed
                    self.printer.write("Updating %s psipred dataframes" % len(seqbuddy))
                    psipred_dfs = rdmcl.update_psipred(alignment, psipred_dfs, "trimal")

                # Prepare all-by-all list
                self.printer.write("Preparing all-by-all data")
                data_len, data = rdmcl.prepare_all_by_all(seqbuddy, psipred_dfs, self.cpus)

                if num_subjobs == 1 and data_len > self.cpus * self.job_size_coff:
                    data_len, data, subjob_num, num_subjobs = self.spawn_subjobs(id_hash, data, psipred_dfs,
                                                                                 gap_open, gap_extend)
                elif subjob_num > 1:
                    data_len, data = self.load_subjob(id_hash, subjob_num, num_subjobs, psipred_dfs)

                # Launch multicore
                self.printer.write("Running all-by-all data (%s comparisons)" % data_len)
                with open(self.data_file, "w") as ofile:
                    ofile.write("seq1,seq2,subsmat,psi")

                br.run_multicore_function(data, rdmcl.mc_score_sequences, quiet=True, max_processes=self.cpus,
                                          func_args=[alignment, gap_open, gap_extend, self.data_file])

                self.printer.write("Processing final results")
                self.process_final_results(id_hash, subjob_num, num_subjobs)

                self.running += time.time() - self.split_time
                self.split_time = time.time()

            except (OSError, FileNotFoundError, br.GuessError, ValueError) as err:
                if num_subjobs == 1:
                        self.terminate("something wrong with primary cluster %s\n%s" % (full_name, err))
                else:
                    with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                        cursor.execute("DELETE FROM processing WHERE hash=?", (full_name,))
                    continue

        # Broken out of while loop, clean up and terminate worker
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)

        self.terminate("deleted check file")

    def idle_workers(self):
        # Calculate how many workers are just sitting around
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            workers = cursor.execute("SELECT COUNT(*) FROM heartbeat WHERE thread_type='worker'").fetchone()[0]

        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            jobs = cursor.execute("SELECT * FROM processing WHERE hash NOT LIKE '%%^_%%' ESCAPE '^'").fetchall()
            subjobs = cursor.execute("SELECT * FROM processing WHERE hash LIKE '%%^_%%' ESCAPE '^'").fetchall()

        jobs = [x[0] for x in jobs]
        subjobs = [x[0] for x in subjobs]
        active_jobs = len(subjobs)
        subjobs = set([x.split("_")[-1] for x in subjobs])
        for job in jobs:
            active_jobs += 1 if job not in subjobs else 0
        return workers - active_jobs if active_jobs < workers else 0  # This can happen if a worker died

    def check_masters(self, idle):
        if time.time() - self.last_heartbeat_from_master > self.max_wait:
            terminate = False
            with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                cursor.execute("SELECT * FROM heartbeat WHERE thread_type='master' "
                               "AND pulse>?", (time.time() - self.max_wait - cursor.lag,))
                masters = cursor.fetchall()
                if not masters:
                    terminate = True
                self.last_heartbeat_from_master = time.time()
            if terminate:
                self.terminate("%s of master inactivity (spent %s%% time idle)" %
                               (br.pretty_time(self.max_wait), idle))
        return

    def clean_dead_threads(self):
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            wait_time = time.time() - self.dead_thread_wait - cursor.lag
            dead_masters = cursor.execute("SELECT thread_id FROM heartbeat WHERE thread_type='master' "
                                          "AND pulse < ?", (wait_time,)).fetchall()
            if dead_masters:
                dead_masters = ", ".join([str(x[0]) for x in dead_masters])
                cursor.execute("DELETE FROM heartbeat WHERE thread_id IN (%s)" % dead_masters)
            else:
                dead_masters = ""

            master_ids = cursor.execute("SELECT thread_id FROM heartbeat WHERE thread_type='master'").fetchall()
            master_ids = "" if not master_ids else ", ".join([str(x[0]) for x in master_ids])

            dead_workers = cursor.execute("SELECT * FROM heartbeat WHERE thread_type='worker' "
                                          "AND pulse < ?", (wait_time,)).fetchall()
            if dead_workers:
                dead_workers = ", ".join([str(x[0]) for x in dead_workers])
                cursor.execute("DELETE FROM heartbeat WHERE thread_id IN (%s)" % dead_workers)

        with helpers.ExclusiveConnect(self.wrkdb_path, max_lock=120) as cursor:
            # Remove any jobs in the 'processing' table where the worker is dead
            if dead_workers:
                cursor.execute("DELETE FROM processing WHERE worker_id IN (%s)" % dead_workers)
            # Add master ids from orphaned entries to the dead masters list
            orphans = cursor.execute("SELECT master_id FROM waiting "
                                     "WHERE master_id NOT IN (%s)" % master_ids).fetchall()
            orphans += cursor.execute("SELECT master_id FROM proc_comp "
                                      "WHERE master_id NOT IN (%s)" % master_ids).fetchall()
            if orphans:
                orphans = ", ".join([str(x[0]) for x in orphans])
                dead_masters += ", %s" % orphans if dead_masters else orphans

            cursor.execute("DELETE FROM waiting WHERE master_id IN (%s)" % dead_masters)
            cursor.execute("DELETE FROM proc_comp WHERE master_id IN (%s)" % dead_masters)
            waiting_hashes = cursor.execute("SELECT hash FROM waiting").fetchall()
            waiting_hashes = [x[0] for x in waiting_hashes]
            dead_hashes = []
            all_hashes = cursor.execute("SELECT hash FROM queue").fetchall()
            all_hashes += cursor.execute("SELECT hash FROM processing").fetchall()
            all_hashes += cursor.execute("SELECT hash FROM complete").fetchall()
            all_hashes = [x[0] for x in all_hashes]
            for next_hash in all_hashes:
                if next_hash.split("_")[-1] not in waiting_hashes:
                    dead_hashes.append(next_hash)
            if dead_hashes:
                dead_hashes = set(dead_hashes)
                for id_hash in dead_hashes:
                    for del_file in ["%s.%s" % (id_hash, x) for x in ["graph", "aln", "seqs"]]:
                        try:
                            os.remove(os.path.join(self.output, del_file))
                        except FileNotFoundError:
                            pass
                        shutil.rmtree(os.path.join(self.output, id_hash), ignore_errors=True)
                dead_hashes = "'%s'" % "', '".join(dead_hashes)
                cursor.execute("DELETE FROM queue WHERE hash IN (%s)" % dead_hashes)
                cursor.execute("DELETE FROM complete WHERE hash IN (%s)" % dead_hashes)
        return

    def fetch_queue_job(self):
        while True:
            with helpers.ExclusiveConnect(self.wrkdb_path, priority=True) as cursor:
                data = cursor.execute('SELECT * FROM queue').fetchone()
                if data:
                    id_hash, psipred_dir, align_m, align_p, trimal, gap_open, gap_extend = data
                    if trimal:
                        trimal = trimal.split()
                        for indx, arg in enumerate(trimal):
                            try:
                                trimal[indx] = float(arg)
                            except ValueError:
                                pass

                    cursor.execute("DELETE FROM queue WHERE hash=?", (id_hash,))
                    if cursor.execute('SELECT worker_id FROM processing WHERE hash=?', (id_hash,)).fetchone():
                        continue

                    cursor.execute("INSERT INTO processing (hash, worker_id)"
                                   " VALUES (?, ?)", (id_hash, self.heartbeat.id,))
                    data = [id_hash, psipred_dir, align_m, align_p, trimal, gap_open, gap_extend]
                break
        return data

    def prepare_psipred_dfs(self, seqbuddy, psipred_dir):
        psipred_dfs = OrderedDict()
        self.printer.write("Preparing %s psipred dataframes" % len(seqbuddy))
        for rec in seqbuddy.records:
            psipred_file = os.path.join(psipred_dir, "%s.ss2" % rec.id)
            if not os.path.isfile(psipred_file):
                raise FileNotFoundError(psipred_file)
            psipred_dfs[rec.id] = rdmcl.read_ss2_file(psipred_file)
        return psipred_dfs

    def process_final_results(self, id_hash, subjob_num, num_subjobs):
        with open(self.data_file, "r") as ifile:
            sim_scores = pd.read_csv(ifile, index_col=False)

        if num_subjobs > 1:
            # If all subjobs are not complete, an empty df is returned
            sim_scores = self.process_subjob(id_hash, sim_scores, subjob_num, num_subjobs)

        if not sim_scores.empty:
            sim_scores = rdmcl.set_final_sim_scores(sim_scores)

            with helpers.ExclusiveConnect(os.path.join(self.output, "write.lock"), max_lock=0):
                # Place these write commands in ExclusiveConnect to ensure a writing lock
                sim_scores.to_csv(os.path.join(self.output, "%s.graph" % id_hash), header=None, index=False)

            with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                # Confirm that the job is still being waited on and wasn't killed before adding to the `complete` table
                waiting = cursor.execute("SELECT master_id FROM waiting WHERE hash=?", (id_hash,)).fetchall()
                processing = cursor.execute("SELECT worker_id FROM processing WHERE hash=?",
                                            (id_hash,)).fetchall()
                if waiting and processing:
                    cursor.execute("INSERT INTO complete (hash) "
                                   "VALUES (?)", (id_hash,))
                elif not waiting:
                    for del_file in ["%s.%s" % (id_hash, x) for x in ["graph", "aln", "seqs"]]:
                        try:
                            os.remove(os.path.join(self.output, del_file))
                        except FileNotFoundError:
                            pass

                cursor.execute("DELETE FROM processing WHERE hash=?", (id_hash,))
                cursor.execute("DELETE FROM complete WHERE hash LIKE '%%_%s'" % id_hash)
            if os.path.isdir(os.path.join(self.output, id_hash)):
                shutil.rmtree(os.path.join(self.output, id_hash))
        return

    def spawn_subjobs(self, id_hash, data, psipred_dfs, gap_open, gap_extend):
        subjob_out_dir = os.path.join(self.output, id_hash)
        os.makedirs(subjob_out_dir, exist_ok=True)
        # Flatten the data jobs list back down
        data = [y for x in data for y in x]
        len_data = len(data)

        for rec_id, df in psipred_dfs.items():
            df.to_csv(os.path.join(subjob_out_dir, "%s.ss2" % rec_id), header=None, index=False, sep=" ")

        # Break it up again into min number of chunks where len(each chunk) < #CPUs * job_size_coff
        num_subjobs = int(rdmcl.ceil(len_data / (self.cpus * self.job_size_coff)))
        job_size = int(rdmcl.ceil(len_data / num_subjobs))
        data = [data[i:i + job_size] for i in range(0, len_data, job_size)]

        for indx, subjob in enumerate(data):
            with open(os.path.join(subjob_out_dir, "%s_of_%s.txt" % (indx + 1, num_subjobs)), "w") as ofile:
                for pair in subjob:
                    ofile.write("%s %s\n" % (pair[0], pair[1]))

        # Push all jobs into the queue except the first job, which is held back for the current node to work on
        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            for indx, subjob in enumerate(data[1:]):
                # NOTE: the 'indx + 2' is necessary to push index to '1' start and account for the job already removed
                cursor.execute("INSERT INTO queue (hash, psi_pred_dir, gap_open, gap_extend) "
                               "VALUES (?, ?, ?, ?)", ("%s_%s_%s" % (indx + 2, num_subjobs, id_hash),
                                                       subjob_out_dir, gap_open, gap_extend,))

            cursor.execute("INSERT INTO processing (hash, worker_id) VALUES (?, ?)",
                           ("1_%s_%s" % (num_subjobs, id_hash), self.heartbeat.id,))

        n = int(rdmcl.ceil(len(data[0]) / self.cpus))
        data = [data[0][i:i + n] for i in range(0, len(data[0]), n)]

        subjob_num = 1
        return len(data[0]), data, subjob_num, num_subjobs

    def load_subjob(self, id_hash, subjob_num, num_subjobs, psipred_dfs):
        subjob_out_dir = os.path.join(self.output, id_hash)
        with open(os.path.join(subjob_out_dir, "%s_of_%s.txt" % (subjob_num, num_subjobs)), "r") as ifile:
            data = ifile.read().strip().split("\n")

        data = [line.split() for line in data]
        data = [(rec1, rec2, psipred_dfs[rec1], psipred_dfs[rec2]) for rec1, rec2 in data]
        data_len = len(data)
        n = int(rdmcl.ceil(len(data) / self.cpus))
        data = [data[i:i + n] for i in range(0, data_len, n)]
        return data_len, data

    def process_subjob(self, id_hash, sim_scores, subjob_num, num_subjobs):
        output = pd.DataFrame()
        subjob_out_dir = os.path.join(self.output, id_hash)
        if os.path.isdir(subjob_out_dir):
            full_id_hash = "%s_%s_%s" % (subjob_num, num_subjobs, id_hash)
            sim_scores.to_csv(os.path.join(subjob_out_dir, "%s_of_%s.sim_df" % (subjob_num, num_subjobs)), index=False)

            with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                # Confirm that the job is still being waited on and wasn't killed before adding to the `complete` table
                waiting = cursor.execute("SELECT master_id FROM waiting WHERE hash=?", (id_hash,)).fetchall()
                processing = cursor.execute("SELECT worker_id FROM processing WHERE hash=? AND worker_id=?",
                                            (full_id_hash, self.heartbeat.id,)).fetchall()
                complete = cursor.execute("SELECT hash FROM complete WHERE hash=?", (full_id_hash,)).fetchall()

                if waiting and processing and not complete:
                    cursor.execute("INSERT INTO complete (hash) "
                                   "VALUES (?)", (full_id_hash,))

                cursor.execute("DELETE FROM processing WHERE hash=?", (full_id_hash,))
                complete_count = cursor.execute("SELECT COUNT(*) FROM complete "
                                                "WHERE hash LIKE '%%_%s'" % id_hash).fetchone()[0]

            if complete_count == num_subjobs:
                output = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi"])
                for indx in range(1, num_subjobs + 1):
                    next_df = os.path.join(subjob_out_dir, "%s_of_%s.sim_df" % (indx, num_subjobs))
                    output = output.append(pd.read_csv(next_df, index_col=False))
        return output

    def terminate(self, message):
        self.printer.write("Terminating Worker_%s because of %s." % (self.heartbeat.id, message))
        self.printer.new_line(1)
        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            cursor.execute("DELETE FROM processing WHERE worker_id=?", (self.heartbeat.id,))
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)
        if os.path.isfile(self.worker_file):
            os.remove(self.worker_file)
        self.heartbeat.end()
        sys.exit()


# Used to patch br.DynamicPrint() when the -log flag is thrown
def _write(self):
    try:
        while True:
            self.out_type.write("\n%s" % self._next_print,)
            self.out_type.flush()
            self._last_print = self._next_print
            yield
    finally:
        pass


def argparse_init():
    import argparse

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="launch_worker", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                     description='''\
\033[1mLaunch Worker\033[m
  Call the union, it's time to put those idle nodes to work.

  Set up a queue and launch a worker node to calculate 
  all-by-all similarity graphs for RD-MCL.

\033[1mUsage\033[m:
  launch_worker [-options]
''')

    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(), metavar="",
                              help="Specify the directory where sqlite databases will be fed by RD-MCL", )
    parser_flags.add_argument("-hr", "--heart_rate", type=int, default=60, metavar="",
                              help="Specify how often the worker should check in (default=60)")
    parser_flags.add_argument("-mw", "--max_wait", action="store", type=int, default=600, metavar="",
                              help="Specify the maximum time a worker will stay alive "
                                   "without seeing a master (default=600)")
    parser_flags.add_argument("-dtw", "--dead_thread_wait", action="store", type=int, default=120, metavar="",
                              help="Specify the maximum time a worker will wait to see a "
                                   "heartbeat before killing a thread (default=120)")
    parser_flags.add_argument("-cpu", "--max_cpus", type=int, action="store", default=cpu_count(), metavar="",
                              help="Specify the maximum number of cores the worker can use (default=%s)" % cpu_count())
    parser_flags.add_argument("-js", "--job_size", type=int, action="store", default=cpu_count(), metavar="",
                              help="Set job size coffactor to adjust… well, job size. (default=300)")
    parser_flags.add_argument("-log", "--log", help="Stream log data one line at a time", action="store_true")
    parser_flags.add_argument("-q", "--quiet", help="Suppress all output", action="store_true")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()

    workdb = os.path.join(in_args.workdb, "work_db.sqlite")
    heartbeatdb = os.path.join(in_args.workdb, "heartbeat_db.sqlite")

    connection = sqlite3.connect(workdb)
    cur = connection.cursor()
    for sql in ['CREATE TABLE queue (hash TEXT PRIMARY KEY, psi_pred_dir TEXT, '
                'align_m TEXT, align_p TEXT, trimal TEXT, gap_open FLOAT, gap_extend FLOAT)',
                'CREATE TABLE processing (hash TEXT PRIMARY KEY, worker_id INTEGER)',
                'CREATE TABLE complete   (hash TEXT PRIMARY KEY)',
                'CREATE TABLE proc_comp   (hash TEXT PRIMARY KEY, master_id INTEGER)',
                'CREATE TABLE waiting (hash TEXT, master_id INTEGER)']:
        try:
            cur.execute(sql)
        except sqlite3.OperationalError:
            pass

    cur.close()
    connection.close()

    connection = sqlite3.connect(heartbeatdb)
    cur = connection.cursor()
    try:
        cur.execute('CREATE TABLE heartbeat (thread_id INTEGER PRIMARY KEY AUTOINCREMENT, '
                    'thread_type TEXT, pulse INTEGER)')
    except sqlite3.OperationalError:
        pass

    cur.close()
    connection.close()

    wrkr = Worker(in_args.workdb, heartrate=in_args.heart_rate, max_wait=in_args.max_wait,
                  dead_thread_wait=in_args.dead_thread_wait, cpus=in_args.max_cpus,
                  job_size_coff=in_args.job_size, log=in_args.log, quiet=in_args.quiet)
    valve = br.SafetyValve(5)
    while True:  # The only way out is through Worker.terminate
        try:
            valve.step("Too many Worker crashes detected.")
            wrkr.start()

        except KeyboardInterrupt:
            wrkr.terminate("KeyboardInterrupt")

        except Exception as err:
            with helpers.ExclusiveConnect(wrkr.wrkdb_path) as cursor:
                cursor.execute("DELETE FROM processing WHERE worker_id=?", (wrkr.heartbeat.id,))

            if "Too many Worker crashes detected" in str(err):
                wrkr.terminate("too many Worker crashes")

            tb = "%s: %s\n\n" % (type(err).__name__, err)
            for _line in traceback.format_tb(sys.exc_info()[2]):
                if os.name == "nt":
                    _line = re.sub('"(?:[A-Za-z]:)*\{0}.*\{0}(.*)?"'.format(os.sep), r'"\1"', _line)
                else:
                    _line = re.sub('"{0}.*{0}(.*)?"'.format(os.sep), r'"\1"', _line)
                tb += _line
            print("\nWorker_%s crashed!\n" % wrkr.heartbeat.id, tb)


if __name__ == '__main__':
    main()
