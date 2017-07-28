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
import time
from random import randint, random
from collections import OrderedDict
from io import StringIO
from copy import copy
# My packages
try:
    import helpers
    import rdmcl
except ImportError:
    from . import helpers
    from . import rdmcl

# Globals
WORKERLOCK = Lock()
CPUS = cpu_count()
printer = br.DynamicPrint()


class Worker(object):
    def __init__(self, wrkdb_path, hbdb_path, heartrate=60, max_wait=120):
        self.wrkdb_path = wrkdb_path
        self.hbdb_path = hbdb_path
        self.heartrate = heartrate
        self.last_heartbeat = self.heartrate + time.time()
        self.max_wait = max_wait
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('worker', %s)" % round(time.time()))
            self.id = cursor.lastrowid
        with open("Worker_%s" % self.id, "w") as ofile:
            ofile.write("To terminate this Worker, simply delete this file.")
        self.cpus = br.cpu_count() - 1
        self.data_file = ".Worker_%s.dat" % self.id
        self.start_time = time.time()
        self.split_time = 0
        self.idle = 1
        self.running = 1

    def start(self):
        self.split_time = time.time()
        self.start_time = time.time()
        printer.write("Starting Worker_%s" % self.id)
        printer.new_line(1)
        seqs, psipred_dir, master_id = ["", "", 0]  # Instantiate some variables

        idle_countdown = 1
        while os.path.isfile("Worker_%s" % self.id):
            idle = round(100 * self.idle / (self.idle + self.running), 2)
            if not idle_countdown:
                printer.write("Idle %s%%" % idle)
                idle_countdown = 5

            if time.time() > self.last_heartbeat:
                # Make sure there are some masters still kicking around
                self.last_heartbeat = self.heartrate + time.time()
                with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                    cursor.execute('UPDATE heartbeat SET pulse=%s WHERE thread_id=%s' % (round(time.time()), self.id))
                    cursor.execute("SELECT * FROM heartbeat WHERE thread_type='master' "
                                   "AND pulse>%s" % (time.time() - self.max_wait))
                    masters = cursor.fetchall()
                    if not masters:
                        printer.write("Terminating Worker_%s due to lack of activity by any master threads.\n"
                                      "Spent %s%% of time idle." % (self.id, idle))
                        printer.new_line(1)
                        cursor.execute("DELETE FROM heartbeat WHERE thread_id=%s" % self.id)
                        break

            # Check for and clean up dead threads and orphaned jobs every hundredth(ish) time through
            rand_check = random()
            if rand_check > 0.99:
                with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                    dead_masters = cursor.execute("SELECT * FROM heartbeat WHERE thread_type='master' "
                                                  "AND pulse < %s" % (time.time() - self.max_wait)).fetchall()
                    dead_workers = cursor.execute("SELECT * FROM heartbeat WHERE thread_type='worker' "
                                                  "AND pulse < %s" % (time.time() - self.max_wait)).fetchall()
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
                            if not cursor.execute("SELECT hash FROM waiting WHERE hash='%s'" % job_hash).fetchall():
                                cursor.execute("DELETE FROM complete WHERE hash='%s'" % job_hash)

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

            with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                cursor.execute('SELECT * FROM queue')
                data = cursor.fetchone()
                if data:
                    id_hash, seqs, psipred_dir, master_id = data
                    cursor.execute("INSERT INTO processing (hash, worker_id, master_id)"
                                   " VALUES ('%s', '%s', %s)" % (id_hash, self.id, master_id))
                    cursor.execute("DELETE FROM queue WHERE hash='%s'" % id_hash)

            self.idle += time.time() - self.split_time
            self.split_time = time.time()
            if not data:
                time.sleep(random() * 3)  # Pause for some part of three seconds
                idle_countdown -= 1
                continue

            idle_countdown = 1
            # Prepare alignment
            seqbuddy = Sb.SeqBuddy(StringIO(seqs))
            if len(seqbuddy) == 1:
                alignment = Alb.AlignBuddy(str(seqbuddy))
            else:
                printer.write("Creating MSA (%s seqs)" % len(seqbuddy))
                alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft",
                                             params="--globalpair --thread -1", quiet=True)

            # Prepare psipred dataframes
            psipred_dfs = OrderedDict()
            breakout = False
            printer.write("Preparing %s psipred dataframes" % len(seqbuddy))
            for rec in alignment.records_iter():
                psipred_file = "%s/%s.ss2" % (psipred_dir, rec.id)
                if not os.path.isfile(psipred_file):
                    printer.write("Terminating Worker_%s because psi file %s not found." % (self.id, psipred_file))
                    printer.new_line(1)
                    with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                        cursor.execute('DELETE FROM heartbeat WHERE thread_id=%s' % self.id)

                    with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                        cursor.execute("DELETE FROM processing WHERE hash='%s'" % id_hash)
                    breakout = True
                    break
                psipred_dfs[rec.id] = rdmcl.read_ss2_file(psipred_file)
            if breakout:
                break

            # Need to specify what columns the PsiPred files map to now that there are gaps.
            for rec in seqbuddy.records:
                ss_file = psipred_dfs[rec.id]
                ss_counter = 0
                for indx, residue in enumerate(rec.seq):
                    if residue != "-":
                        psipred_dfs[rec.id].set_value(ss_counter, "indx", indx)
                        ss_counter += 1
                psipred_dfs[rec.id] = ss_file

            # Scores seem to be improved by removing gaps. ToDo: Need to test this explicitly for the paper
            printer.write("Trimal (%s seqs)" % len(seqbuddy))
            alignment = Alb.trimal(alignment, threshold=0.9)

            # Re-update PsiPred files now that some columns, possibly including non-gap characters, are removed
            printer.write("Updating %s psipred dataframes" % len(seqbuddy))
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

            # Prepare all-by-all list
            printer.write("Preparing all-by-all data")
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
            n = int(rdmcl.ceil(len(data) / self.cpus))
            data = [data[i:i + n] for i in range(0, len(data), n)]
            # Launch multicore
            printer.write("Running all-by-all data (%s comparisons)" % data_len)
            with open(".Worker_%s.dat" % self.id, "w") as ofile:
                ofile.write("seq1,seq2,subsmat,psi")

            br.run_multicore_function(data, score_sequences, quiet=True, max_processes=self.cpus,
                                      func_args=[alignment, -5, 0, ".Worker_%s.dat" % self.id, self.id, self.hbdb_path])

            printer.write("Processing final results")
            with open(".Worker_%s.dat" % self.id, "r") as ifile:
                sim_scores = pd.read_csv(ifile, index_col=False)

            # Set raw score, which is used by Orphan placement
            sim_scores['raw_score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
            # Distribute final scores_components between 0-1.
            sim_scores['psi'] = (sim_scores['psi'] - sim_scores['psi'].min()) / \
                                (sim_scores['psi'].max() - sim_scores['psi'].min())

            sim_scores['subsmat'] = (sim_scores['subsmat'] - sim_scores['subsmat'].min()) / \
                                    (sim_scores['subsmat'].max() - sim_scores['subsmat'].min())

            # ToDo: Experiment testing these magic number weights...
            sim_scores['score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)

            with helpers.ExclusiveConnect(self.wrkdb_path) as cursor:
                cursor.execute("INSERT INTO complete (hash, alignment, graph, worker_id, master_id) "
                               "VALUES ('%s', '%s', '%s', %s, %s)" % (id_hash, str(alignment),
                                                                      sim_scores.to_csv(header=None, index=False),
                                                                      self.id, master_id))

                cursor.execute("DELETE FROM processing WHERE hash='%s'" % id_hash)

            self.running += time.time() - self.split_time
            self.split_time = time.time()

        if os.path.isfile(self.data_file):
            os.remove(self.data_file)

        if os.path.isfile("Worker_%s" % self.id):
            os.remove("Worker_%s" % self.id)
        else:
            printer.write("Terminating Worker_%s because check file was deleted." % self.id)
            printer.new_line(1)
            with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                cursor.execute('DELETE FROM heartbeat WHERE thread_id=%s' % self.id)
        return


def score_sequences(data, func_args):
    results = ["" for _ in data]
    alb_obj, gap_open, gap_extend, output_file, worker_id, hbdb_path = func_args
    for indx, recs in enumerate(data):
        align = Alb.make_copy(alb_obj)
        id1, id2, psi1_df, psi2_df = recs

        # Occasionally update the heartbeat database so masters know there is still life
        if random() > 0.99:
            with helpers.ExclusiveConnect(hbdb_path) as cursor:
                cursor.execute('UPDATE heartbeat SET pulse=%s WHERE thread_id=%s' % (round(time.time()), worker_id))

        # Calculate the best possible scores, and divide by the observed scores
        id_regex = "^%s$|^%s$" % (id1, id2)

        # Alignment comparison
        Alb.pull_records(align, id_regex)

        observed_score = 0
        observed_len = align.lengths()[0]
        seq1_best = 0
        seq1_len = 0
        seq2_best = 0
        seq2_len = 0
        seq1, seq2 = align.records()
        prev_aa1 = "-"
        prev_aa2 = "-"

        for aa_pos in range(observed_len):
            aa1 = seq1.seq[aa_pos]
            aa2 = seq2.seq[aa_pos]

            if aa1 != "-":
                seq1_best += rdmcl.BLOSUM62[aa1, aa1]
                seq1_len += 1
            if aa2 != "-":
                seq2_best += rdmcl.BLOSUM62[aa2, aa2]
                seq2_len += 1
            if aa1 == "-" or aa2 == "-":
                if prev_aa1 == "-" or prev_aa2 == "-":
                    observed_score += gap_extend
                else:
                    observed_score += gap_open
            else:
                observed_score += rdmcl.BLOSUM62[aa1, aa2]
            prev_aa1 = str(aa1)
            prev_aa2 = str(aa2)

        # Calculate average per-residue log-odds ratios for both best possible alignments and observed
        # Note: The best score range is 4 to 11. Possible observed range is -4 to 11.
        observed_score = (observed_score / observed_len)

        seq1_best = (seq1_best / seq1_len)
        seq1_score = observed_score / seq1_best

        seq2_best = (seq2_best / seq2_len)
        seq2_score = observed_score / seq2_best

        subs_mat_score = (seq1_score + seq2_score) / 2

        # PSI PRED comparison
        align = Sb.SeqBuddy(align.records())
        Sb.clean_seq(align)
        num_gaps = 0
        ss_score = 0
        for row1 in psi1_df.itertuples():
            row2 = psi2_df.loc[psi2_df["indx"] == row1.indx]
            if not row2.empty:
                row_score = 0
                row_score += 1 - abs(float(row1.coil_prob) - float(row2.coil_prob))
                row_score += 1 - abs(float(row1.helix_prob) - float(row2.helix_prob))
                row_score += 1 - abs(float(row1.sheet_prob) - float(row2.sheet_prob))
                ss_score += row_score / 3
            else:
                num_gaps += 1
        align_len = len(psi2_df) + num_gaps
        ss_score /= align_len
        results[indx] = "\n%s,%s,%s,%s" % (id1, id2, subs_mat_score, ss_score)

    with WORKERLOCK:
        with open(output_file, "a") as ofile:
            ofile.write("".join(results))
    return


def main():
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

    workdb = os.path.join(in_args.workdb, "work_db.sqlite")
    heartbeatdb = os.path.join(in_args.workdb, "heartbeat_db.sqlite")

    global printer
    if in_args.log:
        def _write(self):
            try:
                while True:
                    self.out_type.write("\n%s" % self._next_print,)
                    self.out_type.flush()
                    self._last_print = self._next_print
                    yield
            finally:
                pass
        br.DynamicPrint._write = _write
        printer = br.DynamicPrint()

    if in_args.quiet:
        printer = br.DynamicPrint(quiet=True)

    connection = sqlite3.connect(workdb)
    cur = connection.cursor()
    for sql in ['CREATE TABLE queue (hash TEXT PRIMARY KEY, seqs TEXT, psi_pred_dir TEXT, master_id INTEGER)',
                'CREATE TABLE processing (hash TEXT PRIMARY KEY, worker_id INTEGER, master_id INTEGER)',
                'CREATE TABLE complete   (hash TEXT PRIMARY KEY, alignment TEXT, graph TEXT, '
                'worker_id INTEGER, master_id INTEGER)',
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

    wrkr = Worker(workdb, heartbeatdb, in_args.heart_rate, in_args.max_wait)
    wrkr.start()

if __name__ == '__main__':
    main()
