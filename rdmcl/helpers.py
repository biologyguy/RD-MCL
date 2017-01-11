#!/usr/bin/env python3
import sqlite3
import sys
import re
import json
import logging
import shutil
from math import log
from time import time
from copy import copy
from hashlib import md5
from multiprocessing import SimpleQueue, Process, Pipe
from buddysuite.buddy_resources import pretty_time, TempDir


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
