#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from .. import helpers
import os
import buddysuite.SeqBuddy
import buddysuite.buddy_resources as br
from hashlib import md5
import sqlite3
from multiprocessing.queues import SimpleQueue
from multiprocessing import Pipe, Process


def test_sqlitebroker_init(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    assert broker.db_file == "%s%sdb.sqlite" % (tmpdir.path, hf.sep)
    assert type(broker.connection) == sqlite3.Connection
    assert type(broker.cursor) == sqlite3.Cursor
    assert type(broker.broker_queue) == SimpleQueue
    assert broker.broker is None


def test_sqlitebroker_create_table(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])
    connect = sqlite3.connect("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    cursor = connect.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    response = cursor.fetchone()
    assert response == ("foo",)
    # Try to create the table again so the method skips through the try block
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])


def test_sqlitebroker_broker_loop(hf, monkeypatch, capsys):
    class MockBrokerLoopGet(object):
        def __init__(self, pipe, modes, sql="SELECT name FROM sqlite_master WHERE type='table'"):
            self.mode = self._get(modes, sql)
            self.sendpipe = pipe

        def _get(self, modes, sql):
            for _mode in modes:
                yield {'mode': _mode, 'sql': sql, 'pipe': self.sendpipe}

        def get(self):
            return next(self.mode)

    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])

    recvpipe, sendpipe = Pipe(False)
    get_dict = {'mode': 'stop', 'sql': "SELECT name FROM sqlite_master WHERE type='table'", 'pipe': sendpipe}
    broker.broker_queue.put(get_dict)
    simple_queue_get = MockBrokerLoopGet(sendpipe, ["sql", "stop"])

    monkeypatch.setattr(SimpleQueue, 'get', simple_queue_get.get)
    broker._broker_loop(broker.broker_queue)

    connect = sqlite3.connect("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    cursor = connect.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    response = cursor.fetchone()
    assert response == ("foo",)

    # Test errors
    simple_queue_get = MockBrokerLoopGet(sendpipe, ["foo_bar"])
    monkeypatch.setattr(SimpleQueue, 'get', simple_queue_get.get)

    with pytest.raises(RuntimeError) as err:
        broker._broker_loop(broker.broker_queue)
    assert "Broker instruction 'foo_bar' not understood." in str(err)

    simple_queue_get = MockBrokerLoopGet(sendpipe, ["sql"], "NONSENSE SQL COMMAND")
    monkeypatch.setattr(SimpleQueue, 'get', simple_queue_get.get)

    with pytest.raises(sqlite3.OperationalError) as err:
        broker._broker_loop(broker.broker_queue)
    assert 'sqlite3.OperationalError: near "NONSENSE": syntax error' in str(err)
    out, err = capsys.readouterr()
    assert "Failed query: NONSENSE SQL COMMAND" in out


def test_sqlitebroker_start_and_stop_broker(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    assert broker.broker is None
    broker.start_broker()
    assert type(broker.broker) == Process
    assert broker.broker.is_alive()

    broker.stop_broker()
    assert not broker.broker.is_alive()


def test_sqlitebroker_query(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])
    with pytest.raises(RuntimeError) as err:
        broker.query("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")
    assert "Broker not running." in str(err)

    broker.start_broker()
    query = broker.query("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")
    assert query == []

    broker.close()
    connect = sqlite3.connect("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    cursor = connect.cursor()
    cursor.execute("SELECT * FROM foo")
    response = cursor.fetchone()
    assert response == (0, 'hello', 25)


def test_sqlitebroker_close(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    assert broker.broker is None
    broker.start_broker()
    assert broker.broker.is_alive()
    broker.close()
    assert not broker.broker.is_alive()


def test_logger(hf):
    tmp = br.TempFile()
    logger = helpers.Logger(tmp.path)
    assert type(logger.logger) == helpers.logging.RootLogger
    assert type(logger.console) == helpers.logging.StreamHandler
    assert logger.logger.level == 20
    assert len(logger.logger.handlers) == 2
    assert type(logger.logger.handlers[1]) == helpers.logging.StreamHandler
    assert logger.console.level == 30

    helpers.logging.info("Some info")
    helpers.logging.warning("Some Warnings")

    logger.move_log("%sfirst.log" % tmp.path)

    with open("%sfirst.log" % tmp.path, "r") as ofile:
        assert ofile.read() == "Some info\nSome Warnings\n"

