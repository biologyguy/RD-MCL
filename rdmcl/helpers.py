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
                if query['mode'] == 'insert':
                    self.cursor.execute("INSERT INTO %s %s VALUES %s" %
                                        (query["table"], query["fields"], query["data"]))
                elif query['mode'] == 'update':
                    self.cursor.execute("UPDATE %s SET %s WHERE %s" % (query["table"], query["data"], query["primary"]))
                elif query['mode'] == 'fetch':
                    hash_id, field, pipe = query[1:4]
                    self.cursor.execute("SELECT (%s) FROM %s WHERE %s" % (query['field'], table, hash_id))
                    response = self.cursor.fetchone()
                    if response is None or len(response) == 0:
                        response = None
                    else:
                        response = response[0]
                    pipe.send(response)
                elif query['mode'] == 'sql':
                    pipe = query['pipe']
                    self.cursor.execute(query['sql'])
                    response = self.cursor.fetchall()
                    #response = [x for x in response]
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

    def push(self, table, ):
        self.broker_queue.put(('push', sql))
        return

    def insert(self, table, data):
        """
        Run the INSERT command on the database
        :param table: Table being added to
        :param data: List of tuples; e.g., [('id', 1234), ('price', 2.45)]
        :type data: list of tuples
        :return:
        """
        # Add single quotes around data strings
        for indx, data_pair in enumerate(data):
            if type(data_pair[1]) == str:
                data[indx] = (data_pair[0], "'%s'" % re.sub("'", "", data_pair[1]))

        fields = "(%s)" % ", ".join([x[0] for x in data])
        data = "(%s)" % ", ".join([x[1] for x in data])
        self.broker_queue.put({'mode': 'insert', 'table': table, 'fields': fields, 'data': data})
        return

    def update(self, table, primary_key, data):
        """
        Run the UPDATE command on the database
        :param table: Table being added to
        :param primary_key: Tuple of primary key field and value, e.g., ('seq_id', 12345)
        :param data: List of tuples; e.g., [('id', 1234), ('price', 2.45)]
        :type data: list of tuples
        :return:
        """
        # Add single quotes around data strings
        for indx, data_pair in enumerate(data):
            if type(data_pair[1]) == str:
                data[indx] = (data_pair[0], "'%s'" % re.sub("'", "", data_pair[1]))

        primary_key = (primary_key[0], "'%s'" % re.sub("'", "", primary_key[1])) \
            if type(primary_key[1]) == str else primary_key

        primary_key = "%s=%s" % (primary_key[0], primary_key[1])
        data = "%s" % ", ".join(["%s=%s" % (x, y) for x, y in data])
        self.broker_queue.put({'mode': 'update', 'table': table,
                               'primary': primary_key, 'data': data})
        return

    def fetch(self, sql):
        recvpipe, sendpipe = Pipe(False)
        self.broker_queue.put({'mode': 'sql', 'sql': sql, 'pipe': sendpipe})
        response = json.loads(recvpipe.recv())
        return response

    def close(self):
        self.stop_broker()
        self.connection.close()


if __name__ == '__main__':
    broker = SQLiteBroker()
    broker.create_table("rdmcl", ["hash TEXT PRIMARY KEY", "seq_ids TEXT"])
    broker.start_broker()
    while True:
        _input = input("New record: ")
        if _input == "quit":
            break
        _input = _input.split(" ")
        #broker.insert("rdmcl", [("hash", _input[0]), ("seq_ids", ",".join(_input[1:]))])
        #broker.update("rdmcl", ("hash", _input[0]), [("seq_ids", ",".join(_input[1:]))])
        print(broker.fetch("SELECT * FROM 'rdmcl' WHERE hash='%s'" % _input[0]))

    broker.close()
