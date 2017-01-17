#!/usr/bin/env python3
import sqlite3
import sys
import re
import json
import logging
import shutil
import pandas as pd
import numpy as np
from math import log, sqrt
from time import time
from copy import copy
from hashlib import md5
from multiprocessing import SimpleQueue, Process, Pipe
from buddysuite.buddy_resources import pretty_time, TempDir, SafetyValve


class SQLiteBroker(object):
    """
    Multithread broker to query a SQLite db
    """
    def __init__(self, db_file="sqlite_db.sqlite"):
        self.db_file = db_file
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
        self.broker_queue = SimpleQueue()
        self.broker = None

    def create_table(self, table_name, fields):
        """
        Make a new table in the database
        :param table_name: What do you want your table named?
        :param fields: field names and any SQL modifiers like type or key commands
        (e.g., ['table_id INT PRIMARY KEY', 'some_data TEXT', 'price INT'])
        :type fields: list
        :return:
        """
        fields = ", ".join(fields)
        try:
            self.cursor.execute("CREATE TABLE %s (%s)" % (table_name, fields))
        except sqlite3.OperationalError:
            pass
        return

    def _broker_loop(self, queue):  # The queue must be passed in explicitly because the process is being spun off
        while True:
            if not queue.empty():
                query = queue.get()
                if query['mode'] == 'sql':
                    pipe = query['pipe']
                    try:
                        self.cursor.execute(query['sql'])
                    except sqlite3.OperationalError as err:
                        print("Failed query: %s" % query['sql'])
                        raise err
                    response = self.cursor.fetchall()
                    pipe.send(json.dumps(response))
                elif query['mode'] == 'stop':
                    break
                else:
                    raise RuntimeError("Broker instruction '%s' not understood." % query['mode'])
                self.connection.commit()

    def start_broker(self):
        if not self.broker:
            self.broker = Process(target=self._broker_loop, args=[self.broker_queue])
            self.broker.daemon = True
            self.broker.start()
        return

    def stop_broker(self):
        self.broker_queue.put({'mode': 'stop'})
        while self.broker.is_alive():
            pass  # Don't move on until the broker is all done doing whatever it might be doing
        return

    def query(self, sql):
        if not self.broker:
            raise RuntimeError("Broker not running. Use the 'start_broker()' method before calling query().")

        recvpipe, sendpipe = Pipe(False)
        self.broker_queue.put({'mode': 'sql', 'sql': sql, 'pipe': sendpipe})
        response = json.loads(recvpipe.recv())
        return response

    def close(self):
        self.stop_broker()
        self.connection.close()
        return


class Logger(object):
    def __init__(self, location):
        self.location = location

        # Set up logging. Use 'info' to write to file only, anything higher will go to both terminal and file.
        logging.basicConfig(filename=self.location, level=logging.INFO, format="")
        self.logger = logging.getLogger()
        self.console = logging.StreamHandler()
        self.console.setLevel(logging.WARNING)
        self.logger.addHandler(self.console)

    def move_log(self, location):
        shutil.move(self.location, location)
        logging.basicConfig(filename=location, level=logging.INFO, format="")
        self.location = location
        return


class Timer(object):
    def __init__(self):
        self.start = round(time())
        self.split_time = round(time())

    def split(self, prefix="", postfix=""):
        split = round(time()) - self.split_time
        self.split_time = round(time())
        return "%s%s%s" % (prefix, pretty_time(split), postfix)

    def total_elapsed(self, prefix="", postfix=""):
        return "%s%s%s" % (prefix, pretty_time(round(time()) - self.start), postfix)


def md5_hash(in_str):
    in_str = str(in_str).encode("utf-8")
    return md5(in_str).hexdigest()


def make_full_mat(subsmat):
    for key in copy(subsmat):
        try:
            # don't over-write the reverse keys if they are already initialized
            subsmat[(key[1], key[0])]
        except KeyError:
            subsmat[(key[1], key[0])] = subsmat[key]
    return subsmat


def bit_score(raw_score):
    """
    :param raw_score: Sum of values pulled from BLOSUM62 substitution matrix
    :return:
    """
    # These values were empirically determined for BLOSUM62 by Altschul
    bit_k_value = 0.035
    bit_lambda = 0.252
    bits = ((bit_lambda * raw_score) - (log(bit_k_value))) / log(2)
    return bits


def markov_clustering(data, inflation, edge_sim_threshold=0.):
    def compare(matrix1, matrix2):
        dif_mat = matrix1 - matrix2
        dif_mat **= 2
        return dif_mat.sum().sum()

    def normalize(matrix):
        for indx in range(len(matrix)):
            column = matrix[:, indx]
            column_sum = column.sum()
            if column_sum != 0:
                matrix[:, indx] = column / column_sum
        return matrix

    def mcl_step(matrix):
        # Expand
        matrix = matrix.dot(matrix)
        # Inflate
        matrix = matrix ** inflation
        # Re-normalize
        matrix = normalize(matrix)
        return matrix

    def apply_edge_sim_threshold(threshold, matrix):
        matrix[matrix <= threshold] = 1e-308
        print(matrix)
        return matrix

    def df_to_transition_matrix(graph_df):
        size = (sqrt(8 * len(graph_df) + 1) + 1) / 2
        if not size.is_integer():
            raise ValueError("The provided dataframe is not a symmetric graph")
        score_sum = graph_df.score.sum()
        graph_df['transition_prob'] = graph_df.score / score_sum
        size = int(size)
        tran_mat = np.zeros([size, size])
        name_order = sorted(list(set(graph_df.seq1.tolist() + graph_df.seq2.tolist())))
        for indx, row in graph_df.iterrows():
            seq1 = name_order.index(row.seq1)
            seq2 = name_order.index(row.seq2)
            tran_mat[seq1][seq2] = row.score
            tran_mat[seq2][seq1] = row.score
        tran_mat = apply_edge_sim_threshold(edge_sim_threshold, tran_mat)
        tran_mat = normalize(tran_mat)
        return tran_mat

    valve = SafetyValve(global_reps=10)
    transition_matrix = df_to_transition_matrix(data)
    sub_states = [pd.DataFrame(transition_matrix)]
    counter = 1
    while True:
        valve.step()
        transition_matrix = mcl_step(transition_matrix)
        sub_states.append(pd.DataFrame(transition_matrix))
        if compare(sub_states[-2], sub_states[-1]) == 0:
            break
        counter += 1
    return sub_states


if __name__ == '__main__':
    sample_df = pd.read_csv("../workshop/mcl/complete_all_by_all.scores", sep="\t", header=None, index_col=False)
    sample_df.columns = ["seq1", "seq2", "score"]
    mcl_steps = markov_clustering(sample_df, 8, 0.6)
    print(mcl_steps[0])
    print(mcl_steps[-1])
