#!/usr/bin/env python3
import sqlite3
from multiprocessing import SimpleQueue, Process, Pipe
import sys
import re
import json


class SQLiteBroker(object):
    """
    Simple multithread broker to insert/update records in a db or to retrieve records
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

    def _broker_loop(self, queue):
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
        recvpipe, sendpipe = Pipe(False)
        self.broker_queue.put({'mode': 'sql', 'sql': sql, 'pipe': sendpipe})
        response = json.loads(recvpipe.recv())
        return response

    def close(self):
        self.stop_broker()
        self.connection.close()
