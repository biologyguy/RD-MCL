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
from copy import copy
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
CPUS = cpu_count()
JOB_SIZE_COFACTOR = 300


class Worker(object):
    def __init__(self, location, heartrate=60, max_wait=120, log=False, quiet=False):
        self.wrkdb_path = os.path.join(location, "work_db.sqlite")
        self.hbdb_path = os.path.join(location, "heartbeat_db.sqlite")
        self.output = os.path.join(location, ".worker_output")

        self.masterclear_path = os.path.split(self.wrkdb_path)[0]
        self.masterclear_path = os.path.join(self.masterclear_path, "MasterClear")
        self.heartrate = heartrate
        self.heartbeat = rdmcl.HeartBeat(self.hbdb_path, self.heartrate, thread_type="worker")
        self.max_wait = max_wait
        self.cpus = br.cpu_count() - 1
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

        with open("Worker_%s" % self.heartbeat.id, "w") as ofile:
            ofile.write("To terminate this Worker, simply delete this file.")
        self.data_file = ".Worker_%s.dat" % self.heartbeat.id

        helpers.dummy_func()

        self.last_heartbeat_from_master = time.time()
        self.printer.write("Starting Worker_%s" % self.heartbeat.id)
        self.printer.new_line(1)

        idle_countdown = 1
        while os.path.isfile("Worker_%s" % self.heartbeat.id):
            idle = round(100 * self.idle / (self.idle + self.running), 2)
            if not idle_countdown:
                self.printer.write("Idle %s%%" % idle)
                idle_countdown = 5

            # Make sure there are some masters still kicking around
            self.check_masters(idle)

            # Check for MasterClear override (file with a number in it, specifying how many seconds is considered dead)
            max_wait = self.masterclear()

            # Check for and clean up dead threads and orphaned jobs every hundredth(ish) time through
            rand_check = random()
            if rand_check > 0.99 or max_wait != self.max_wait:
                self.clean_dead_threads(max_wait)

            # Fetch a job from the queue
            data = self.fetch_queue_job()
            if data:
                full_name, psipred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend = data
                subjob_num, num_subjobs, id_hash = [1, 1, full_name] if len(full_name.split("_")) == 1 \
                    else full_name.split("_")
                subjob_num = int(subjob_num)
                num_subjobs = int(num_subjobs)
                self.printer.write("Running %s" % full_name)
            else:
                time.sleep(random() * 3)  # Pause for some part of three seconds
                idle_countdown -= 1
                self.idle += time.time() - self.split_time
                self.split_time = time.time()
                continue

            idle_countdown = 1
            seqbuddy = Sb.SeqBuddy("%s/%s.seqs" % (self.output, id_hash), in_format="fasta")

            # Prepare alignment
            if len(seqbuddy) == 1:
                alignment = Alb.AlignBuddy(str(seqbuddy), in_format="fasta")
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
            if not psipred_dfs:
                break

            if num_subjobs == 1:  # This is starting a full job from scratch, not a sub-job
                # Need to specify what columns the PsiPred files map to now that there are gaps.
                psipred_dfs = self.update_psipred(alignment, psipred_dfs, "msa")

                # TrimAl
                self.printer.write("Trimal (%s seqs)" % len(seqbuddy))
                alignment = self.trimal(seqbuddy, trimal, alignment)

                with helpers.ExclusiveConnect(os.path.join(self.output, "write.lock"), max_lock=0):
                    # Place these write commands in ExclusiveConnect to ensure a writing lock
                    if not os.path.isfile(os.path.join(self.output, "%s.aln" % id_hash)):
                        alignment.write(os.path.join(self.output, "%s.aln" % id_hash), out_format="fasta")

                # Re-update PsiPred files now that some columns, possibly including non-gap characters, are removed
                self.printer.write("Updating %s psipred dataframes" % len(seqbuddy))
                psipred_dfs = self.update_psipred(alignment, psipred_dfs, "trimal")

            # Prepare all-by-all list
            self.printer.write("Preparing all-by-all data")
            data_len, data = self.prepare_all_by_all(seqbuddy, psipred_dfs)

            if num_subjobs == 1 and data_len > self.cpus * JOB_SIZE_COFACTOR:
                data_len, data, subjob_num, num_subjobs = self.spawn_subjobs(id_hash, data, psipred_dfs, master_id,
                                                                             gap_open, gap_extend)
            elif subjob_num > 1:
                data_len, data = self.load_subjob(id_hash, subjob_num, num_subjobs, psipred_dfs)

            # Launch multicore
            self.printer.write("Running all-by-all data (%s comparisons)" % data_len)
            with open(".Worker_%s.dat" % self.heartbeat.id, "w") as ofile:
                ofile.write("seq1,seq2,subsmat,psi")

            br.run_multicore_function(data, score_sequences, quiet=True, max_processes=self.cpus,
                                      func_args=[alignment, gap_open, gap_extend, ".Worker_%s.dat" % self.heartbeat.id])

            self.printer.write("Processing final results")
            self.process_final_results(id_hash, alignment, master_id, subjob_num, num_subjobs)

            self.running += time.time() - self.split_time
            self.split_time = time.time()

        # Broken out of while loop, clean up and terminate worker
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)

        if os.path.isfile("Worker_%s" % self.heartbeat.id):
            os.remove("Worker_%s" % self.heartbeat.id)
        else:
            self.terminate("deleted check file")
        return

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

    def masterclear(self):
        max_wait = int(self.max_wait)
        with WORKERLOCK:
            # Check MasterClear signal (file in working dir with # of second specified for master heartbeat)
            if os.path.isfile(self.masterclear_path):
                with open(self.masterclear_path, "r") as ifile:
                    try:
                        max_wait = int(ifile.read())
                        self.printer.write("MasterClear signal %s" % max_wait)
                        self.printer.new_line(1)
                    except ValueError as err:
                        self.printer.write(str(err))
                        self.printer.new_line(1)
                os.remove(self.masterclear_path)
            else:
                pass
        return max_wait

    def clean_dead_threads(self, max_wait):
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            wait_time = time.time() - max_wait - cursor.lag
            dead_masters = cursor.execute("SELECT * FROM heartbeat WHERE thread_type='master' "
                                          "AND pulse < ?", (wait_time,)).fetchall()
            dead_workers = cursor.execute("SELECT * FROM heartbeat WHERE thread_type='worker' "
                                          "AND pulse < ?", (wait_time,)).fetchall()
            if dead_masters:
                dead_masters = [str(x[0]) for x in dead_masters]
                dead_masters = ", ".join(dead_masters)
                cursor.execute("DELETE FROM heartbeat WHERE thread_id IN (%s)" % dead_masters)
            if dead_workers:
                dead_workers = [str(x[0]) for x in dead_workers]
                dead_workers = ", ".join(dead_workers)
                cursor.execute("DELETE FROM heartbeat WHERE thread_id IN (%s)" % dead_workers)

            master_ids = cursor.execute("SELECT thread_id FROM heartbeat WHERE thread_type='master'").fetchall()

        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            if dead_masters:
                cursor.execute("DELETE FROM queue WHERE master_id IN (%s)" % dead_masters)
                cursor.execute("DELETE FROM waiting WHERE master_id IN (%s)" % dead_masters)

                complete_jobs = cursor.execute("SELECT hash FROM complete").fetchall()
                for job_hash in complete_jobs:
                    if not cursor.execute("SELECT hash FROM waiting WHERE hash=?", (job_hash[0],)).fetchall():
                        cursor.execute("DELETE FROM complete WHERE hash=?", (job_hash[0],))

            if dead_workers:
                cursor.execute("DELETE FROM processing WHERE worker_id IN (%s)" % dead_workers)

            if master_ids:
                master_ids = ", ".join([str(x[0]) for x in master_ids])
                orphaned_jobs = cursor.execute("SELECT hash FROM complete "
                                               "WHERE master_id NOT IN (%s)" % master_ids).fetchall()
                if orphaned_jobs:
                    orphaned_job_hashes = "'%s'" % "', '".join([x[0] for x in orphaned_jobs])

                    waiting = cursor.execute("SELECT hash FROM waiting "
                                             "WHERE hash IN (%s)" % orphaned_job_hashes).fetchall()
                    orphaned_job_hashes = "'%s'" % "', '".join([x[0] for x in waiting])
                    orphaned_jobs = cursor.execute("SELECT hash FROM complete "
                                                   "WHERE hash NOT IN (%s)" % orphaned_job_hashes).fetchall()
                    orphaned_job_hashes = "'%s'" % "', '".join([x[0] for x in orphaned_jobs])
                    cursor.execute("DELETE FROM complete WHERE hash IN (%s)" % orphaned_job_hashes)
        return

    def fetch_queue_job(self):
        with helpers.ExclusiveConnect(self.wrkdb_path, priority=True) as cursor:
            data = cursor.execute('SELECT * FROM queue').fetchone()
            if data:
                id_hash, psipred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend = data
                if trimal:
                    trimal = trimal.split()
                    for indx, arg in enumerate(trimal):
                        try:
                            trimal[indx] = float(arg)
                        except ValueError:
                            pass

                cursor.execute("INSERT INTO processing (hash, worker_id, master_id)"
                               " VALUES (?, ?, ?)", (id_hash, self.heartbeat.id, master_id,))
                cursor.execute("DELETE FROM queue WHERE hash=?", (id_hash,))
                data = [id_hash, psipred_dir, master_id, align_m, align_p, trimal, gap_open, gap_extend]
        return data

    def prepare_psipred_dfs(self, seqbuddy, psipred_dir):
        psipred_dfs = OrderedDict()
        self.printer.write("Preparing %s psipred dataframes" % len(seqbuddy))
        for rec in seqbuddy.records:
            psipred_file = os.path.join(psipred_dir, "%s.ss2" % rec.id)
            if not os.path.isfile(psipred_file):
                self.terminate("missing psi ss2 file (%s)" % psipred_file)
            psipred_dfs[rec.id] = rdmcl.read_ss2_file(psipred_file)
        return psipred_dfs

    @staticmethod
    def update_psipred(alignment, psipred_dfs, mode):
        if mode == "msa":
            for rec in alignment.records_iter():
                ss_file = psipred_dfs[rec.id]
                ss_counter = 0
                for indx, residue in enumerate(rec.seq):
                    if residue != "-":
                        psipred_dfs[rec.id].set_value(ss_counter, "indx", indx)
                        ss_counter += 1
                psipred_dfs[rec.id] = ss_file

        elif mode == "trimal":
            for rec in alignment.records_iter():
                # Instantiate list of max possible size
                new_psi_pred = [0 for _ in range(len(psipred_dfs[rec.id].index))]
                indx = 0
                for row in psipred_dfs[rec.id].itertuples():
                    if alignment.alignments[0].position_map[int(row[1])][1]:
                        new_psi_pred[indx] = list(row)[1:]
                        indx += 1
                new_psi_pred = new_psi_pred[:indx]
                psipred_dfs[rec.id] = pd.DataFrame(new_psi_pred, columns=["indx", "aa", "ss", "coil_prob",
                                                                          "helix_prob", "sheet_prob"])
        else:
            raise ValueError("Unrecognized mode '%s': select from ['msa', 'trimal']" % mode)
        return psipred_dfs

    @staticmethod
    def trimal(seqbuddy, trimal, alignment):
        # Scores seem to be improved by removing gaps. ToDo: Need to test this explicitly for the paper
        # Only remove columns up to a 50% reduction in average seq length and only if all sequences are retained
        ave_seq_length = Sb.ave_seq_length(seqbuddy)
        for threshold in trimal:
            align_copy = Alb.trimal(Alb.make_copy(alignment), threshold=threshold)
            cleaned_seqs = Sb.clean_seq(Sb.SeqBuddy(str(align_copy)))
            cleaned_seqs = Sb.delete_small(cleaned_seqs, 1)
            # Structured this way for unit test purposes
            if len(alignment.records()) != len(cleaned_seqs):
                continue
            elif Sb.ave_seq_length(cleaned_seqs) / ave_seq_length < 0.5:
                continue
            else:
                alignment = align_copy
                break
        return alignment

    def prepare_all_by_all(self, seqbuddy, psipred_dfs):
        ids1 = [rec.id for rec in seqbuddy.records]
        ids2 = copy(ids1)
        data = [0 for _ in range(int((len(ids1)**2 - len(ids1)) / 2))]
        indx = 0
        for rec1 in ids1:
            del ids2[ids2.index(rec1)]
            for rec2 in ids2:
                data[indx] = (rec1, rec2, psipred_dfs[rec1], psipred_dfs[rec2])
                indx += 1

        data_len = len(data)
        n = int(rdmcl.ceil(data_len / self.cpus))
        data = [data[i:i + n] for i in range(0, data_len, n)]
        return data_len, data

    def process_final_results(self, id_hash, alignment, master_id, subjob_num, num_subjobs):
        with open(".Worker_%s.dat" % self.heartbeat.id, "r") as ifile:
            sim_scores = pd.read_csv(ifile, index_col=False)

        if num_subjobs > 1:
            sim_scores = self.process_subjob(id_hash, sim_scores, subjob_num, num_subjobs, master_id)

        if not sim_scores.empty:
            # Set raw score, which is used by Orphan placement
            sim_scores['raw_score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
            # Distribute final scores_components between 0-1.
            sim_scores['psi'] = (sim_scores['psi'] - sim_scores['psi'].min()) / \
                                (sim_scores['psi'].max() - sim_scores['psi'].min())

            sim_scores['subsmat'] = (sim_scores['subsmat'] - sim_scores['subsmat'].min()) / \
                                    (sim_scores['subsmat'].max() - sim_scores['subsmat'].min())

            # ToDo: Experiment testing these magic number weights...
            sim_scores['score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)

            with helpers.ExclusiveConnect(os.path.join(self.output, "write.lock"), max_lock=0):
                # Place these write commands in ExclusiveConnect to ensure a writing lock
                if not os.path.isfile(os.path.join(self.output, "%s.graph" % id_hash)):
                    sim_scores.to_csv(os.path.join(self.output, "%s.graph" % id_hash), header=None, index=False)

            with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                # Confirm that the job is still being waited on before adding to the `complete` table
                waiting = cursor.execute("SELECT master_id FROM waiting WHERE hash=?", (id_hash,)).fetchall()

                if waiting:
                    cursor.execute("INSERT INTO complete (hash, worker_id, master_id) "
                                   "VALUES (?, ?, ?)", (id_hash, self.heartbeat.id, master_id,))
                else:
                    os.remove(os.path.join(self.output, "%s.graph" % id_hash))
                    os.remove(os.path.join(self.output, "%s.aln" % id_hash))
                    os.remove(os.path.join(self.output, "%s.seqs" % id_hash))

                cursor.execute("DELETE FROM processing WHERE hash=?", (id_hash,))
                cursor.execute("DELETE FROM processing WHERE hash LIKE '%%_%s'" % id_hash)
                cursor.execute("DELETE FROM complete WHERE hash LIKE '%%_%s'" % id_hash)
            if os.path.isdir(os.path.join(self.output, id_hash)):
                shutil.rmtree(os.path.join(self.output, id_hash))
        return

    def spawn_subjobs(self, id_hash, data, psipred_dfs, master_id, gap_open, gap_extend):
        out_dir = os.path.join(self.output, id_hash)
        os.makedirs(out_dir, exist_ok=True)
        # Flatten the data jobs list back down
        data = [y for x in data for y in x]
        len_data = len(data)

        for rec_id, df in psipred_dfs.items():
            df.to_csv(os.path.join(out_dir, "%s.ss2" % rec_id), header=None, index=False, sep=" ")

        # Break it up again into min number of chunks where len(each chunk) < #CPUs * JOB_SIZE_COFACTOR
        num_subjobs = int(rdmcl.ceil(len_data / (self.cpus * JOB_SIZE_COFACTOR)))
        job_size = int(rdmcl.ceil(len_data / num_subjobs))
        data = [data[i:i + job_size] for i in range(0, len_data, job_size)]

        for indx, subjob in enumerate(data):
            with open(os.path.join(out_dir, "%s_of_%s.txt" % (indx + 1, num_subjobs)), "w") as ofile:
                for pair in subjob:
                    ofile.write("%s %s\n" % (pair[0], pair[1]))

        # Push all jobs into the queue except the first job, which is held back for the current node to work on
        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            for indx, subjob in enumerate(data[1:]):
                # NOTE: the 'indx + 2' is necessary to push index to '1' start and account for the job already removed
                cursor.execute("INSERT INTO queue (hash, psi_pred_dir, master_id, gap_open, gap_extend) "
                               "VALUES (?, ?, ?, ?, ?)", ("%s_%s_%s" % (indx + 2, num_subjobs, id_hash),
                                                          out_dir, master_id, gap_open, gap_extend,))

            cursor.execute("INSERT INTO processing (hash, worker_id, master_id) VALUES (?, ?, ?)",
                           ("1_%s_%s" % (num_subjobs, id_hash), self.heartbeat.id, master_id,))

        n = int(rdmcl.ceil(len(data[0]) / self.cpus))
        data = [data[0][i:i + n] for i in range(0, len(data[0]), n)]

        subjob_num = 1
        return len(data[0]), data, subjob_num, num_subjobs

    def load_subjob(self, id_hash, subjob_num, num_subjobs, psipred_dfs):
        out_dir = os.path.join(self.output, id_hash)
        with open(os.path.join(out_dir, "%s_of_%s.txt" % (subjob_num, num_subjobs)), "r") as ifile:
            data = ifile.read().strip().split("\n")

        data = [line.split() for line in data]
        data = [(rec1, rec2, psipred_dfs[rec1], psipred_dfs[rec2]) for rec1, rec2 in data]
        data_len = len(data)
        n = int(rdmcl.ceil(len(data) / self.cpus))
        data = [data[i:i + n] for i in range(0, data_len, n)]
        return data_len, data

    def process_subjob(self, id_hash, sim_scores, subjob_num, num_subjobs, master_id):
        out_dir = os.path.join(self.output, id_hash)
        full_id_hash = "%s_%s_%s" % (subjob_num, num_subjobs, id_hash)
        sim_scores.to_csv(os.path.join(out_dir, "%s_of_%s.sim_df" % (subjob_num, num_subjobs)), index=False)

        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            cursor.execute("INSERT INTO complete (hash, worker_id, master_id) "
                           "VALUES (?, ?, ?)", (full_id_hash, self.heartbeat.id, master_id,))
            cursor.execute("DELETE FROM processing WHERE hash=?", (full_id_hash,))
            complete_count = cursor.execute("SELECT COUNT(*) FROM complete "
                                            "WHERE hash LIKE '%%_%s'" % id_hash).fetchone()[0]

        sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi"])
        if complete_count == num_subjobs:
            for indx in range(1, num_subjobs + 1):
                next_df = os.path.join(out_dir, "%s_of_%s.sim_df" % (indx, num_subjobs))
                sim_scores = sim_scores.append(pd.read_csv(next_df, index_col=False))
        else:
            sim_scores = pd.DataFrame()

        return sim_scores

    def terminate(self, message):
        self.printer.write("Terminating Worker_%s because of %s." % (self.heartbeat.id, message))
        self.printer.new_line(1)
        with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
            cursor.execute("DELETE FROM processing WHERE worker_id=?", (self.heartbeat.id,))
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)
        if os.path.isfile("Worker_%s" % self.heartbeat.id):
            os.remove("Worker_%s" % self.heartbeat.id)
        self.heartbeat.end()
        sys.exit()


def score_sequences(data, func_args):
    # ##################################################################### #
    # Calculate the best possible scores, and divide by the observed scores #
    # ##################################################################### #

    results = ["" for _ in data]
    alb_obj, gap_open, gap_extend, output_file = func_args
    for indx, recs in enumerate(data):
        alb_obj_copy = Alb.make_copy(alb_obj)
        id1, id2, psi1_df, psi2_df = recs

        id_regex = "^%s$|^%s$" % (id1, id2)

        # Alignment comparison
        alb_obj_copy = Alb.pull_records(alb_obj_copy, id_regex)

        subs_mat_score = rdmcl.compare_pairwise_alignment(alb_obj_copy, gap_open, gap_extend)

        # PSI PRED comparison
        ss_score = rdmcl.compare_psi_pred(psi1_df, psi2_df)
        results[indx] = "\n%s,%s,%s,%s" % (id1, id2, subs_mat_score, ss_score)

    with WORKERLOCK:
        with open(output_file, "a") as ofile:
            ofile.write("".join(results))
    return


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

    parser = argparse.ArgumentParser(prog="launch_worker", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(),
                        help="Specify the directory where sqlite databases will be fed by RD-MCL", )
    parser.add_argument("-hr", "--heart_rate", type=int, default=60,
                        help="Specify how often the worker should check in")
    parser.add_argument("-mw", "--max_wait", action="store", type=int, default=120,
                        help="Specify the maximum time a worker will stay alive without seeing a master")
    parser.add_argument("-log", "--log", help="Stream log data one line at a time", action="store_true")
    parser.add_argument("-q", "--quiet", help="Suppress all output", action="store_true")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()

    workdb = os.path.join(in_args.workdb, "work_db.sqlite")
    heartbeatdb = os.path.join(in_args.workdb, "heartbeat_db.sqlite")
    worker_output = os.path.join(in_args.workdb, ".worker_output")

    connection = sqlite3.connect(workdb)
    cur = connection.cursor()
    for sql in ['CREATE TABLE queue (hash TEXT PRIMARY KEY, psi_pred_dir TEXT, master_id INTEGER, '
                'align_m TEXT, align_p TEXT, trimal TEXT, gap_open FLOAT, gap_extend FLOAT)',
                'CREATE TABLE processing (hash TEXT PRIMARY KEY, worker_id INTEGER, master_id INTEGER)',
                'CREATE TABLE complete   (hash TEXT PRIMARY KEY, worker_id INTEGER, master_id INTEGER)',
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

    os.makedirs(worker_output, exist_ok=True)

    wrkr = Worker(in_args.workdb, heartrate=in_args.heart_rate, max_wait=in_args.max_wait,
                  log=in_args.log, quiet=in_args.quiet)
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
