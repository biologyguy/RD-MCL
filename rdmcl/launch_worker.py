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
    def __init__(self, wrkdb_path, hbdb_path, heartbeat_wait=120):
        self.wrkdb_path = wrkdb_path
        self.hbdb_path = hbdb_path
        self.heartbeat_wait = heartbeat_wait
        self.last_heartbeat = self.heartbeat_wait + time.time()
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

        while os.path.isfile("Worker_%s" % self.id):
            printer.write("Idle %s%%" % round(100 * self.idle / (self.idle + self.running), 2))
            with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                if time.time() > self.last_heartbeat:
                    cursor.execute("SELECT * FROM heartbeat WHERE thread_type='master'")
                    masters = cursor.fetchall()
                    breakout = True
                    for master in masters:
                        if time.time() < master[2] + self.heartbeat_wait:
                            breakout = False
                            break
                    if breakout:
                        printer.write("Terminating Worker_%s due to lack of activity by any master threads." % self.id)
                        printer.new_line(1)
                        break
                cursor.execute('UPDATE heartbeat SET pulse=%s WHERE thread_id=%s' % (round(time.time()), self.id))

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
                time.sleep(randint(1, 100) / 100)  # Pause for some part of one second
                continue

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
                    print("Terminating Worker_%s because psi file %s not found." % (self.id, psipred_file))
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
                    data[indx] = (rec1, rec2, alignment, psipred_dfs[rec1], psipred_dfs[rec2],
                                  -5, 0, ".Worker_%s.dat" % self.id, self.id, self.hbdb_path)
                    indx += 1

            # Launch multicore
            printer.write("Running all-by-all data (%s comparisons)" % len(data))
            with open(".Worker_%s.dat" % self.id, "w") as ofile:
                ofile.write("seq1,seq2,subsmat,psi")

            br.run_multicore_function(data, score_sequences, quiet=True, max_processes=self.cpus)

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
            print("Terminating Worker_%s because check file was deleted." % self.id)
            with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                cursor.execute('DELETE FROM heartbeat WHERE thread_id=%s' % self.id)
        return


def score_sequences(args):
    id1, id2, align, psi1_df, psi2_df, gap_open, gap_extend, output_file, worker_id, hbdb_path = args

    # Occasionally update the heartbeat database so masters know there is still life
    if random() < 2 / CPUS:
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
    with WORKERLOCK:
        with open(output_file, "a") as ofile:
            ofile.write("\n%s,%s,%s,%s" % (id1, id2, subs_mat_score, ss_score))
    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="launch_worker", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(),
                        help="Specify the directory where sqlite databases will be fed by RD-MCL", )
    parser.add_argument("-mw", "--max_wait", help="", action="store", type=int, default=120)

    in_args = parser.parse_args()

    workdb = os.path.join(in_args.workdb, "work_db.sqlite")
    heartbeatdb = os.path.join(in_args.workdb, "heartbeat_db.sqlite")

    connection = sqlite3.connect(workdb)
    cur = connection.cursor()
    for sql in ['CREATE TABLE queue (hash TEXT PRIMARY KEY, seqs TEXT, psi_pred_dir TEXT, master_id INTEGER)',
                'CREATE TABLE processing (hash TEXT PRIMARY KEY, worker_id INTEGER, master_id INTEGER)',
                'CREATE TABLE complete   (hash TEXT PRIMARY KEY, alignment TEXT, graph TEXT, '
                'worker_id INTEGER, master_id INTEGER)']:
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

    wrkr = Worker(workdb, heartbeatdb, in_args.max_wait)
    wrkr.start()
