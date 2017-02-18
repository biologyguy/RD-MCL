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
from time import time, sleep
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
                        if "database is locked" in str(err):
                            # Wait a few seconds and try one more time, it might get through.
                            sleep(5)
                            self.cursor.execute(query['sql'])
                        else:
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


class MarkovClustering(object):
    def __init__(self, data, inflation, edge_sim_threshold=0.):
        self.dataframe = data
        self.inflation = inflation
        self.edge_sim_threshold = edge_sim_threshold
        self.name_order = sorted(list(set(self.dataframe.seq1.tolist() + self.dataframe.seq2.tolist())))
        self.trans_matrix = self._df_to_transition_matrix()
        self.sub_state_dfs = [pd.DataFrame(self.trans_matrix)]
        self.clusters = []

    @staticmethod
    def compare(df1, df2):
        dif = df1 - df2
        dif **= 2
        return dif.sum().sum()

    @staticmethod
    def normalize(np_matrix):
        for indx in range(len(np_matrix)):
            column = np_matrix[:, indx]
            column_sum = column.sum()
            if column_sum != 0:
                np_matrix[:, indx] = column / column_sum
        return np_matrix

    def _df_to_transition_matrix(self):
        size = (sqrt(8 * len(self.dataframe) + 1) + 1) / 2
        if not size.is_integer():
            raise ValueError("The provided dataframe is not a symmetric graph")
        size = int(size)
        tran_mat = np.zeros([size, size])
        for indx, row in self.dataframe.iterrows():
            seq1 = self.name_order.index(row.seq1)
            seq2 = self.name_order.index(row.seq2)
            tran_mat[seq1][seq2] = row.score
            tran_mat[seq2][seq1] = row.score
        tran_mat[tran_mat <= self.edge_sim_threshold] = 0

        # This is a 'centering' step that is used by the original MCL algorithm
        # ToDo: Compare the outcomes of not centering and of skewing the center to higher values
        for i in range(len(tran_mat)):
            tran_mat[i][i] = max(tran_mat[i])

        tran_mat = self.normalize(tran_mat)
        return tran_mat

    def mcl_step(self):
        # Expand
        # ToDo: There is an issue here, with simulated data and the next command dies quietly
        self.trans_matrix = self.trans_matrix.dot(self.trans_matrix)
        # Inflate
        self.trans_matrix = self.trans_matrix ** self.inflation
        # Re-normalize
        self.trans_matrix = self.normalize(self.trans_matrix)
        return

    def run(self):
        valve = SafetyValve(global_reps=1000)
        while True:
            try:
                valve.step()
            except RuntimeError:
                self.clusters = [self.name_order]
                return
            self.mcl_step()
            self.sub_state_dfs.append(pd.DataFrame(self.trans_matrix))
            if self.compare(self.sub_state_dfs[-2], self.sub_state_dfs[-1]) == 0:
                break

        next_cluster = []
        not_clustered = list(self.name_order)
        for i, row in self.sub_state_dfs[-1].iterrows():
            for j, cell in row.items():
                if cell != 0 and self.name_order[j] in not_clustered:
                    next_cluster.append(self.name_order[j])
                    del not_clustered[not_clustered.index(self.name_order[j])]
            if next_cluster:
                self.clusters.append(next_cluster)
                next_cluster = []
        self.clusters += [[x] for x in not_clustered]
        self.clusters.sort(key=len)
        self.clusters.reverse()
        return

    def write(self, ofile="clusters.mcl"):
        with open(ofile, "w") as _ofile:
            for cluster in self.clusters:
                _ofile.write("%s\n" % "\t".join(cluster))
