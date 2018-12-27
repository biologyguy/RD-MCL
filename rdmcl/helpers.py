#!/usr/bin/env python3
import sqlite3
import json
import logging
import shutil
import os
import re
import pandas as pd
import numpy as np
from scipy import stats
from math import log, sqrt, ceil
from time import time, sleep
from copy import copy
from hashlib import md5
from collections import OrderedDict
from random import random
from multiprocessing import SimpleQueue, Process, Pipe
from subprocess import PIPE, check_output, CalledProcessError

from buddysuite import buddy_resources as br
import signal

RED = "\033[91m"
GREEN = "\033[92m"
DEF_FONT = "\033[39m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
END = '\033[0m'

# Set global precision levels
np.set_printoptions(precision=12)
pd.set_option("display.precision", 12)

contributor_list = [br.Contributor("Stephen", "Bond", commits=520, github="https://github.com/biologyguy"),
                    br.Contributor("Karl", "Keat", commits=30, github="https://github.com/KarlKeat")]

SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))

try:
    git_commit = check_output(['git', '--git-dir={0}{1}..{1}.git'.format(SCRIPT_PATH, os.sep), 'rev-parse',
                               '--short', 'HEAD'], stderr=PIPE).decode().strip()
    git_commit = " (git %s)" % git_commit if git_commit else ""
except CalledProcessError:
    git_commit = ""

VERSION = br.Version("", 1, "2b" + git_commit, contributor_list, {"year": 2018, "month": 1, "day": 3})


class AttrWrapper(object):
    def __init__(self, wrapped):
        self._wrapped = wrapped

    def __getattr__(self, n):
        return getattr(self._wrapped, n)


class ExclusiveConnect(object):
    def __init__(self, db_path, log_message=None, priority=False, log_path="ExclusiveConnect.log", max_lock=60):
        self.db_path = db_path
        self.log_message = log_message
        self.log_output = []
        self.start_time = time()
        self.loop_counter = 0
        self.log_path = log_path
        self.priority = 0.5 if not priority else 1000
        self.max_lock = max_lock

    def raise_timeout(self, *args):
        raise EnvironmentError("ExclusiveConnect Lock held for over %s seconds" % self.max_lock)

    def __enter__(self):
        # Note that there is a pseudo-priority counter
        while True:
            try:
                self.connection = sqlite3.connect(self.db_path, isolation_level=None, timeout=self.priority)
                self.connection.execute('BEGIN EXCLUSIVE')
                if self.log_message:
                    self.log_output.append(round(time() - self.start_time, 4))
                    self.log_output.append(self.loop_counter)
                break
            except sqlite3.OperationalError as err:
                if "database is locked" in str(err):
                    self.loop_counter += 1
                    sleep(1 / self.priority)
                    self.priority += 0.1
                    continue
                else:
                    raise err

        cursor = AttrWrapper(self.connection.cursor())
        cursor.lag = time() - self.start_time
        signal.signal(signal.SIGALRM, self.raise_timeout)
        signal.alarm(self.max_lock)
        return cursor

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.connection.commit()
        self.connection.close()
        if self.log_message:
            self.log_output.append(round(time() - self.start_time, 4))
            self.log_output.append(round(self.log_output[2] - self.log_output[0], 4))
            self.log_output.append(self.log_message)
            log_output = "%s\n" % "\t".join([str(x) for x in self.log_output])
            with open(self.log_path, "a") as ofile:
                ofile.write(log_output)
        signal.alarm(0)


class SQLiteBroker(object):
    """
    Multithread broker to query a SQLite db
    """
    def __init__(self, db_file="sqlite_db.sqlite", lock_wait_time=120):
        self.db_file = db_file
        self.connection = sqlite3.connect(self.db_file)
        self.broker_cursor = self.connection.cursor()
        self.broker_queue = SimpleQueue()
        self.broker = None
        self.lock_wait_time = lock_wait_time
        # ToDo: Set up a process pool to limit number of query threads

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
            self.broker_cursor.execute("CREATE TABLE %s (%s)" % (table_name, fields))
        except sqlite3.OperationalError:
            pass
        return

    def _broker_loop(self, queue):  # The queue must be passed in explicitly because the process is being spun off
        while True:
            if not queue.empty():
                query = queue.get()
                if query['mode'] == 'sql':
                    pipe = query['pipe']
                    locked_counter = 0
                    while True:
                        try:
                            dummy_func()
                            self.broker_cursor.execute(query['sql'], query['values'])
                            self.connection.commit()
                        except sqlite3.OperationalError as err:
                            if "database is locked" in str(err):
                                # Wait for database to become free
                                if locked_counter > self.lock_wait_time * 5:
                                    print("Failed query: %s" % query['sql'])
                                    raise err
                                locked_counter += 1
                                sleep(.2)
                                continue
                            else:
                                print("Failed query: %s" % query['sql'])
                                raise err
                        break
                    response = self.broker_cursor.fetchall()
                    pipe.send(json.dumps(response))
                elif query['mode'] == 'stop':
                    break
                else:
                    raise RuntimeError("Broker instruction '%s' not understood." % query['mode'])

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

    def query(self, sql, values=None, errors=True):
        """
        :param sql: SQL string
        :param values: If question marks are used in SQL command, pass in replacement values as tuple
        :param errors: Suppress raised errors by passing in False
        :return: 
        """
        if not self.broker:
            raise RuntimeError("Broker not running. Use the 'start_broker()' method before calling query().")

        values = () if not values else values
        recvpipe, sendpipe = Pipe(False)
        valve = br.SafetyValve(150)
        while True:
            valve.step("To many threads being called, tried for 5 minutes but couldn't find an open thread.")
            try:
                dummy_func()
                self.broker_queue.put({'mode': 'sql', 'sql': sql, 'values': values, 'pipe': sendpipe})
                break
            except (sqlite3.Error, sqlite3.OperationalError, sqlite3.IntegrityError, sqlite3.DatabaseError) as err:
                if errors:
                    raise err
            except RuntimeError as err:
                if "can't start new thread" in str(err):
                    sleep(2)
                else:
                    raise err
        response = json.loads(recvpipe.recv())
        return response

    def iterator(self, sql):  # Note that this does not run through the broker
        temp_cursor = self.connection.cursor()
        query_result = temp_cursor.execute(sql)
        while True:
            fetched = query_result.fetchone()
            if not fetched:
                break
            else:
                yield fetched

    def close(self):
        self.stop_broker()
        self.connection.close()
        return


class Logger(object):
    def __init__(self, location):
        self.location = location
        open(location, "w").close()

        # Set up logging. Use 'info' to write to file only, anything higher will go to both terminal and file.
        logging.basicConfig(filename=self.location, level=logging.INFO, format="")
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)
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
        return "%s%s%s" % (prefix, br.pretty_time(split), postfix)

    def total_elapsed(self, prefix="", postfix=""):
        return "%s%s%s" % (prefix, br.pretty_time(round(time()) - self.start), postfix)


class KellysColors(object):
    # https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
    def __init__(self):
        self.kelly_colors = OrderedDict(deep_yellowish_brown=(89, 51, 21),
                                        strong_reddish_brown=(127, 24, 13),
                                        strong_purplish_red=(179, 40, 81),
                                        strong_purplish_pink=(246, 118, 142),
                                        vivid_red=(193, 0, 32),
                                        vivid_reddish_orange=(241, 58, 19),
                                        vivid_orange=(255, 104, 0),
                                        strong_yellowish_pink=(255, 122, 92),
                                        vivid_orange_yellow=(255, 142, 0),
                                        vivid_yellow=(255, 179, 0),
                                        vivid_greenish_yellow=(244, 200, 0),
                                        grayish_yellow=(206, 162, 98),
                                        vivid_yellowish_green=(147, 170, 0),
                                        vivid_green=(0, 125, 52),
                                        dark_olive_green=(35, 44, 22),
                                        very_light_blue=(166, 189, 215),
                                        strong_blue=(0, 83, 138),
                                        strong_violet=(83, 55, 122),
                                        strong_purple=(128, 62, 117),
                                        medium_gray=(129, 112, 102))

    def color_iter(self):
        degree = 0
        while True:
            for color, rgb in self.kelly_colors.items():
                rgb = purturb_rgb(rgb, degree)
                yield rgb
            degree += 10


def purturb_rgb(rgb, degree=10):
    new_rgb = []
    for code in rgb:
        degree = round(random() * degree)
        degree = degree * -1 if random() < 0.5 else degree
        while degree != 0:
            code += degree
            if code > 255:
                degree = 255 - code
                code = 255
            elif code < 0:
                degree = abs(code)
                code = 0
            else:
                degree = 0
        new_rgb.append(code)
    return new_rgb[0], new_rgb[1], new_rgb[2]


def mean(series):
    return np.around(np.mean(series), 12)


def std(series):
    return np.around(np.std(series), 12)


def dummy_func(*args, **kwargs):
    """
    This can be placed in code for unit test monkey patching
    :param args: arguments
    :param kwargs: key-word arguments
    :return:
    """
    return args, kwargs


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


def create_truncnorm(mu, sigma, lower=0, upper=1):
    sigma = sigma if sigma > 0.001 else 0.001  # This prevents unrealistically small differences and DivBy0 errors
    dist = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    return dist


def prepare_clusters(ifile,  hierarchy=False):
    with open(ifile, "r") as ifile:
        output = ifile.readlines()
    if output[-1] == "\n":
        del output[-1]

    if hierarchy:
        for indx, line in enumerate(output):
            group_name = re.search("(^.*?)\s", line).group(1)
            line = re.sub(r"(^.*?)\s-*[0-9]+\.[0-9]*\s+", "", line)
            line = line.split()
            output[indx] = (group_name, line)
        output = OrderedDict(output)
    else:
        for indx, line in enumerate(output):
            line = re.sub(r"(^.*?)\s-*[0-9]+\.[0-9]*\s+", "", line)
            line = line.split()
            output[indx] = line
    return output


def chunk_list(l, num_chunks):
    """
    Break up a list into a list of lists
    :param l: Input list
    :param num_chunks: How many lists should the list be chunked into
    :return:
    """
    num_chunks = int(num_chunks)
    if num_chunks < 1 or not l:
        raise AttributeError("Input list must have items in it and num_chunks must be a positive integer")

    size = int(ceil(len(l) / num_chunks))
    num_long = len(l) % num_chunks
    num_long = num_long if num_long != 0 else num_chunks
    chunks = [l[i:i + size] for i in range(0, num_long * size, size)]
    if size != 1:
        chunks += [l[i:i + size - 1] for i in range(num_long * size, len(l), size - 1)]
    return chunks


# ToDo: implement Regularized MCL to take into account flows of neighbors (Expansion step is M*M_G, instead of M*M)
# https://www.youtube.com/watch?v=574z9nisRuE around 12:00
class MarkovClustering(object):
    def __init__(self, data, inflation, edge_sim_threshold=0.):
        """
        Run MCL on a set if data
        :param data: Pandas dataframe of form {"seq1", "seq2", "score"}
        :param inflation: Inflation value
        :param edge_sim_threshold: Make the graph sparse by removing any edges below this threshold
        """
        self.dataframe = data.copy()
        self.dataframe.seq1 = self.dataframe.seq1.astype('category')
        self.dataframe.seq2 = self.dataframe.seq2.astype('category')
        self.dataframe.score = self.dataframe.score.astype('float32')
        self.inflation = inflation
        self.edge_sim_threshold = edge_sim_threshold
        self.name_order = {}
        self.name_order_indx = {}
        for indx, seq_id in enumerate(sorted(list(set(self.dataframe.seq1.tolist() + self.dataframe.seq2.tolist())))):
            self.name_order[seq_id] = indx
            self.name_order_indx[indx] = seq_id

        self.trans_matrix = self._df_to_transition_matrix()
        self.clusters = []

    @staticmethod
    def compare(np1, np2):
        dif = np1 - np2
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
        tran_mat = np.zeros([size, size], dtype="float32")

        for indx1, seq1, seq2, score in self.dataframe[["seq1", "seq2", "score"]].itertuples():
            seq1_indx = self.name_order[seq1]
            seq2_indx = self.name_order[seq2]
            tran_mat[seq1_indx][seq2_indx] = score
            tran_mat[seq2_indx][seq1_indx] = score

        tran_mat[tran_mat <= self.edge_sim_threshold] = 0

        # This is a 'centering' step that is used by the original MCL algorithm
        # ToDo: Compare the outcomes of not centering and of skewing the center to higher values
        for i in range(len(tran_mat)):
            tran_mat[i][i] = max(tran_mat[i])

        tran_mat = self.normalize(tran_mat)
        return tran_mat

    def finalize_transition_matrix(self):
        """
        The transition matrix needs to be reduced to 0s and 1s, with the 1s being read from each row to signify a
        cluster.
        Each column sums to 1.0, and only one sequence per column can be set to the member of that particular cluster,
        so transpose the transition matrix so columns become rows, then find the max score in the row. If the score is
        1.0, then you can just move on to the next row, otherwise, record the max value and the column index to compare
        against later rows, and set the remaining cells to 0. When another row has the same max value other than 1.0 in
        a particular column, then that value is also changed to 1, for inclusion in the group.
        The transition matrix is modified in place, and transposed back to normal before the method returns.
        :return: None
        """
        self.trans_matrix = self.trans_matrix.T
        finished_columns = {}
        for row_indx, row in enumerate(self.trans_matrix):
            cells = sorted([(col_indx, value) for col_indx, value in enumerate(row)], key=lambda x: (x[1], -x[0]),
                           reverse=True)
            if cells[0][1] != 1:
                if cells[0][0] in finished_columns:
                    self.trans_matrix[row_indx][cells[0][0]] = 0. if cells[0][1] < finished_columns[cells[0][0]] else 1.
                else:
                    self.trans_matrix[row_indx][cells[0][0]] = 1.
                    finished_columns[cells[0][0]] = cells[0][1]

                for col_indx, value in cells[1:]:
                    self.trans_matrix[row_indx][col_indx] = 0.
                    if value == 0.:
                        break
        self.trans_matrix = self.trans_matrix.T
        return

    def mcl_step(self):
        # Expand
        self.trans_matrix = np.matmul(self.trans_matrix, self.trans_matrix)
        # Inflate
        self.trans_matrix = self.trans_matrix ** self.inflation
        # Re-normalize
        self.trans_matrix = self.normalize(self.trans_matrix)
        return

    def run(self):
        valve = br.SafetyValve(global_reps=1000)
        last_substate = self.trans_matrix.copy()
        while True:
            try:
                valve.step()
            except RuntimeError:  # No convergence after 1000 MCL steps
                break
            self.mcl_step()
            if self.compare(last_substate, self.trans_matrix) == 0:
                break
            else:
                last_substate = self.trans_matrix.copy()

        if len(np.where(np.logical_and(self.trans_matrix > 0., self.trans_matrix < 1.))[0]):
            self.finalize_transition_matrix()

        next_cluster = []
        not_clustered = list(self.name_order)
        for row in self.trans_matrix:
            for i, cell in enumerate(row):
                if cell != 0 and self.name_order_indx[i] in not_clustered:
                    next_cluster.append(self.name_order_indx[i])
                    del not_clustered[not_clustered.index(self.name_order_indx[i])]
            if next_cluster:
                self.clusters.append(next_cluster)
                next_cluster = []
        self.clusters += [[x] for x in not_clustered]
        self.clusters = sorted(self.clusters, key=lambda x: len(x))
        self.clusters.reverse()
        return

    def write(self, ofile="clusters.mcl"):
        with open(ofile, "w") as _ofile:
            for cluster in self.clusters:
                _ofile.write("%s\n" % "\t".join(cluster))
