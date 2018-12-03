#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from rdmcl import helpers
import os
import buddysuite.buddy_resources as br
import sqlite3
import pandas as pd
import numpy as np
import time
from multiprocessing.queues import SimpleQueue
from multiprocessing import Pipe, Process
from Bio.SubsMat import SeqMat, MatrixInfo
from io import StringIO
from types import SimpleNamespace

pd.set_option('expand_frame_repr', False)


def test_attrwraper():
    test_obj = helpers.AttrWrapper(object())
    test_obj.lag = 3
    assert test_obj.lag == 3


def test_exclusiveconnect(monkeypatch):
    tmpdir = br.TempDir()

    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
    cursor = connect.cursor()
    cursor.execute("CREATE TABLE foo (id INT PRIMARY KEY, some_data TEXT, numbers INT)")
    connect.commit()

    # Basic instantiation
    excl_con = helpers.ExclusiveConnect(os.path.join(tmpdir.path, "db.sqlite"))
    assert excl_con.db_path == os.path.join(tmpdir.path, "db.sqlite")
    assert excl_con.log_message is None
    assert excl_con.log_output == []
    assert excl_con.start_time < time.time()
    assert excl_con.loop_counter == 0
    assert excl_con.log_path == "ExclusiveConnect.log"
    assert excl_con.priority == 0.5
    assert excl_con.max_lock == 60

    # Timeout method
    excl_con.max_lock = 120
    with pytest.raises(EnvironmentError) as err:
        excl_con.raise_timeout()
    assert "ExclusiveConnect Lock held for over 120 seconds" in str(err)

    # Non-logged activity
    with helpers.ExclusiveConnect(os.path.join(tmpdir.path, "db.sqlite")) as cursor:
        cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")

    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
    cursor = connect.cursor()
    cursor.execute("SELECT some_data FROM foo WHERE id=0")
    response = cursor.fetchone()
    assert response == ("hello",)

    # Logged activity
    log_file = tmpdir.subfile("log_file")
    with helpers.ExclusiveConnect(os.path.join(tmpdir.path, "db.sqlite"),
                                  log_message="Testing logging", log_path=log_file) as cursor:
        cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (1, 'bonjour', 50)")

    with open(log_file, "r") as ifile:
        assert "Testing logging" in ifile.read()

    # Locked database with sleep
    class SQLiteError(object):
        def __init__(self):
            self.next_err = self.error_loop()

        @staticmethod
        def error_loop():
            _errors = [sqlite3.OperationalError("database is locked"), sqlite3.OperationalError("Not expected")]
            for _err in _errors:
                yield _err

        def raise_error(self, *args, **kwargs):
            raise next(self.next_err)

    errors = SQLiteError()
    monkeypatch.setattr(sqlite3, "connect", errors.raise_error)
    with pytest.raises(sqlite3.OperationalError) as err:
        with helpers.ExclusiveConnect(os.path.join(tmpdir.path, "db.sqlite")) as cursor:
            cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (2, 'hola', 75)")
    assert "Not expected" in str(err)


def test_sqlitebroker_init():
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    assert broker.db_file == os.path.join(tmpdir.path, "db.sqlite")
    assert type(broker.connection) == sqlite3.Connection
    assert type(broker.broker_cursor) == sqlite3.Cursor
    assert type(broker.broker_queue) == SimpleQueue
    assert broker.broker is None


def test_sqlitebroker_create_table():
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])
    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
    cursor = connect.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    response = cursor.fetchone()
    assert response == ("foo",)
    # Try to create the table again so the method skips through the try block
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])


def test_sqlitebroker_broker_loop(monkeypatch, capsys):
    class MockBrokerLoopGet(object):
        def __init__(self, pipe, modes, sql="SELECT name FROM sqlite_master WHERE type='table'"):
            self.mode = self._get(modes, sql)
            self.sendpipe = pipe

        def _get(self, modes, sql):
            for _mode in modes:
                yield {'mode': _mode, 'sql': sql, 'pipe': self.sendpipe, "values": ()}

        def get(self):
            return next(self.mode)

    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])

    recvpipe, sendpipe = Pipe(False)
    get_dict = {'mode': 'stop', 'sql': "SELECT name FROM sqlite_master WHERE type='table'", 'pipe': sendpipe}
    broker.broker_queue.put(get_dict)
    simple_queue_get = MockBrokerLoopGet(sendpipe, ["sql", "stop"])

    monkeypatch.setattr(SimpleQueue, 'get', simple_queue_get.get)
    broker._broker_loop(broker.broker_queue)

    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
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

    simple_queue_get = MockBrokerLoopGet(sendpipe, ["sql"], "")
    monkeypatch.setattr(SimpleQueue, 'get', simple_queue_get.get)

    def raise_error():
        raise sqlite3.OperationalError("database is locked")
    monkeypatch.setattr(helpers, 'dummy_func', raise_error)
    broker.lock_wait_time = 0.1
    with pytest.raises(sqlite3.OperationalError) as err:
        broker._broker_loop(broker.broker_queue)
    assert 'database is locked' in str(err)
    out, err = capsys.readouterr()
    assert out == 'Failed query: \n'


def test_sqlitebroker_start_and_stop_broker():
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    assert broker.broker is None
    broker.start_broker()
    assert type(broker.broker) == Process
    assert broker.broker.is_alive()

    broker.stop_broker()
    assert not broker.broker.is_alive()


def test_sqlitebroker_query(monkeypatch):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    broker.create_table("foo", ['id INT PRIMARY KEY', 'some_data TEXT', 'numbers INT'])
    with pytest.raises(RuntimeError) as err:
        broker.query("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")
    assert "Broker not running." in str(err)

    broker.start_broker()
    query = broker.query("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")
    assert query == []

    broker.close()
    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
    cursor = connect.cursor()
    cursor.execute("SELECT * FROM foo")
    response = cursor.fetchone()
    assert response == (0, 'hello', 25)

    def raise_error(*_, **__):
        raise sqlite3.OperationalError("sqlite error")

    monkeypatch.setattr(SimpleQueue, "put", raise_error)
    with pytest.raises(sqlite3.OperationalError) as err:
        broker.query("NONSENSE QUERY")

    assert "sqlite error" in str(err)

    def raise_error(*_, **__):
        raise RuntimeError("can't start new thread")

    monkeypatch.setattr(SimpleQueue, "put", raise_error)
    monkeypatch.setattr(helpers, "sleep", raise_error)
    with pytest.raises(RuntimeError) as err:
        broker.query("NONSENSE QUERY")
    assert "can't start new thread" in str(err)

    def raise_error(*_, **__):
        raise RuntimeError("some other runtime error")

    monkeypatch.setattr(SimpleQueue, "put", raise_error)
    monkeypatch.setattr(helpers, "sleep", raise_error)
    with pytest.raises(RuntimeError) as err:
        broker.query("NONSENSE QUERY")
    assert "some other runtime error" in str(err)


def test_sqlitebroker_iterator():
    tmpdir = br.TempDir()

    connect = sqlite3.connect(os.path.join(tmpdir.path, "db.sqlite"))
    cursor = connect.cursor()
    cursor.execute("CREATE TABLE foo (id INT PRIMARY KEY, some_data TEXT, numbers INT)")
    cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (0, 'hello', 25)")
    cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (1, 'bonjour', 50)")
    cursor.execute("INSERT INTO foo (id, some_data, numbers) VALUES (2, 'hola', 75)")
    connect.commit()
    connect.close()

    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    results = broker.iterator("SELECT * FROM foo")
    assert next(results) == (0, 'hello', 25)
    assert next(results) == (1, 'bonjour', 50)
    assert next(results) == (2, 'hola', 75)
    with pytest.raises(StopIteration):
        next(results)


def test_sqlitebroker_close():
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    assert broker.broker is None
    broker.start_broker()
    assert broker.broker.is_alive()
    broker.close()
    assert not broker.broker.is_alive()


def test_logger():
    tmp = br.TempFile()
    logger = helpers.Logger(tmp.path)

    assert type(logger.logger) == helpers.logging.RootLogger
    assert type(logger.console) == helpers.logging.StreamHandler
    assert logger.logger.level == 20
    handlers = [type(handler) for handler in logger.logger.handlers]
    assert len(logger.logger.handlers) == 2, print(handlers)
    assert type(logger.logger.handlers[1]) == helpers.logging.StreamHandler, print(handlers)
    assert logger.console.level == 30

    logger.logger.log(helpers.logging.WARNING, "Some info")
    helpers.logging.warning("Some Warnings")

    logger.move_log("%sfirst.log" % tmp.path)
    with open("%sfirst.log" % tmp.path, "r") as ofile:
        assert ofile.read() == "Some info\nSome Warnings\n"


def test_timer(monkeypatch):
    monkeypatch.setattr(helpers, "time", lambda *_: 1)
    timer = helpers.Timer()
    monkeypatch.setattr(helpers, "time", lambda *_: 2)
    assert timer.split(prefix="start_", postfix="_end") == 'start_1 sec_end'
    monkeypatch.setattr(helpers, "time", lambda *_: 3)
    assert timer.split(prefix="start_", postfix="_end") == 'start_1 sec_end'
    monkeypatch.setattr(helpers, "time", lambda *_: 4)
    assert timer.total_elapsed(prefix="start_", postfix="_end") == 'start_3 sec_end'


def test_mean(hf):
    data = hf.get_data("cteno_sim_scores")
    assert helpers.mean(data.score) == 0.40629959990800002


def test_std(hf):
    data = hf.get_data("cteno_sim_scores")
    assert helpers.std(data.score) == 0.18440242260100001


def test_md5_hash():
    assert helpers.md5_hash("Hello") == "8b1a9953c4611296a827abf8c47804d7"


def test_make_full_mat():
    blosum62 = helpers.make_full_mat(SeqMat(MatrixInfo.blosum62))
    assert blosum62["A", "B"] == -2
    assert blosum62["B", "A"] == -2


def test_bit_score():
    assert helpers.bit_score(100) == 41.192416298119


def test_create_truncnorm():
    dist = helpers.create_truncnorm(0.5, 0.2)
    assert round(dist.cdf(0.6), 12) == 0.693870199384


def test_chunk_list():
    assert helpers.chunk_list(list(range(10)), 4) == [[0, 1, 2], [3, 4, 5], [6, 7], [8, 9]]
    assert helpers.chunk_list(list(range(10)), 10) == [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]
    assert helpers.chunk_list(list(range(10)), 100) == [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]
    assert helpers.chunk_list(list(range(10)), 1) == [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]

    with pytest.raises(AttributeError) as err:
        helpers.chunk_list(list(range(10)), 0)
    assert "Input list must have items in it and num_chunks must be a positive integer" in str(err)

    with pytest.raises(AttributeError) as err:
        helpers.chunk_list([], 10)
    assert "Input list must have items in it and num_chunks must be a positive integer" in str(err)


def test_markov_clustering_init():
    data = """\
Bab\tCfu\t1
Bab\tOma\t1
Bab\tMle\t0
Cfu\tMle\t0
Cfu\tOma\t1
Oma\tMle\t0"""
    sample_df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                            dtype={"seq1": "category", "seq2": "category", "score": "float32"})
    mcl = helpers.MarkovClustering(sample_df, 2, 0.6)
    assert str(mcl.dataframe) == str(sample_df)
    assert mcl.inflation == 2
    assert mcl.edge_sim_threshold == 0.6
    assert sorted(list(mcl.name_order)) == ['Bab', "Cfu", "Mle", "Oma"]
    assert sorted(list(mcl.name_order_indx)) == [0, 1, 2, 3]
    expected = np.matrix([[0.333333343267, 0.333333343267, 0., 0.333333343267],
                          [0.333333343267, 0.333333343267, 0., 0.333333343267],
                          [0.,         0.,         0., 0.],
                          [0.333333343267, 0.333333343267, 0., 0.333333343267]],
                         dtype=str(mcl.trans_matrix.dtype))
    assert np.array_equal(mcl.trans_matrix, expected), print(mcl.trans_matrix)
    assert mcl.clusters == []


def test_markov_clustering_compare():
    data = """\
Bab\tCfu\t1
Bab\tOma\t1
Bab\tMle\t0
Cfu\tMle\t0
Cfu\tOma\t1
Oma\tMle\t0"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    df1 = pd.DataFrame(mcl.trans_matrix)
    df2 = df1.copy()
    assert helpers.MarkovClustering.compare(df1, df2) == 0

    df2[0][0] = 1
    df2[1][3] = 1.5
    compare = helpers.MarkovClustering.compare(df1, df2)
    assert round(float(compare), 1) == 1.8


def test_markov_clustering_normalize():
    matrix = np.matrix([[0., 1., 0., 1.],
                        [1., 0., 0., 1.],
                        [0., 0., 0., 0.],
                        [1., 1., 0., 0.]])
    normalized = helpers.MarkovClustering.normalize(matrix)
    expected = np.matrix([[0., 0.5, 0., 0.5],
                          [0.5, 0., 0., 0.5],
                          [0., 0., 0., 0.],
                          [0.5, 0.5, 0., 0.]])
    assert np.array_equal(normalized, expected), print(normalized)


def test_markov_clustering_df_to_transition_matrix():
    data = """\
Bab\tCfu\t1
Bab\tOma\t1
Bab\tMle\t0
Cfu\tMle\t0
Cfu\tOma\t1
Oma\tMle\t0"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    expected = np.array([[0.333333343267, 0.333333343267, 0., 0.333333343267],
                         [0.333333343267, 0.333333343267, 0., 0.333333343267],
                         [0., 0., 0., 0.],
                         [0.333333343267, 0.333333343267, 0., 0.333333343267]],
                        dtype=str(mcl.trans_matrix.dtype))
    assert np.array_equal(mcl._df_to_transition_matrix(), expected), print(mcl._df_to_transition_matrix())

    mcl.dataframe = mcl.dataframe.ix[1:, :]
    with pytest.raises(ValueError) as err:
        mcl._df_to_transition_matrix()
    assert "The provided dataframe is not a symmetric graph" in str(err)


def test_markov_clustering_finalize_transition_matrix():
    mtrx = np.array([[1., 0., 1., 0.],
                     [0., 0.3, 0., 0.],
                     [0., 0.3, 0., 1.],
                     [0., 0.4, 0., 0.]])

    mcl_obj = SimpleNamespace(trans_matrix=mtrx)
    helpers.MarkovClustering.finalize_transition_matrix(mcl_obj)
    expected = np.array([[1., 0., 1., 0.],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 1.],
                         [0., 1., 0., 0.]])
    assert np.array_equal(mcl_obj.trans_matrix, expected), print(mcl_obj.trans_matrix)

    mtrx = np.array([[1., 0., 1., 0.],
                     [0., 0.2, 0., 0.],
                     [0., 0.5, 0., 1.],
                     [0., 0.3, 0., 0.]])
    mcl_obj = SimpleNamespace(trans_matrix=mtrx)
    helpers.MarkovClustering.finalize_transition_matrix(mcl_obj)
    expected = np.array([[1., 0., 1., 0.],
                         [0., 0., 0., 0.],
                         [0., 1., 0., 1.],
                         [0., 0., 0., 0.]])
    assert np.array_equal(mcl_obj.trans_matrix, expected), print(mcl_obj.trans_matrix)

    mtrx = np.array([[0.5, 0.5, 0., 0.],
                     [0.5, 0.5, 0., 0.],
                     [0., 0., 0.5, 0.5],
                     [0., 0., 0.5, 0.5]])
    mcl_obj = SimpleNamespace(trans_matrix=mtrx)
    helpers.MarkovClustering.finalize_transition_matrix(mcl_obj)
    expected = np.array([[1., 1., 0., 0.],
                         [0., 0., 0., 0.],
                         [0., 0., 1., 1.],
                         [0., 0., 0., 0.]])
    assert np.array_equal(mcl_obj.trans_matrix, expected), print(mcl_obj.trans_matrix)


def test_markov_clustering_mcl_step():
    data = """\
Bab\tCfu\t1
Bab\tOma\t1
Bab\tMle\t0
Cfu\tMle\t1
Cfu\tOma\t1
Oma\tMle\t0"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    mcl.mcl_step()
    expected = np.array([[0.325268805027, 0.197712421417, 0.050000000745, 0.325268805027],
                         [0.325268805027, 0.472222208977, 0.449999988079, 0.325268805027],
                         [0.024193543941, 0.132352933288, 0.449999988079, 0.024193543941],
                         [0.325268805027, 0.197712421417, 0.050000000745, 0.325268805027]],
                        dtype=str(mcl.trans_matrix.dtype))
    assert np.array_equal(mcl.trans_matrix, expected), print(mcl.trans_matrix)

    mcl.mcl_step()
    expected = np.array([[0.256081040, 0.179641350, 0.06652327, 0.256081040],
                         [0.471649080, 0.581160800, 0.64254290, 0.471649080],
                         [0.016188763, 0.059556555, 0.22441064, 0.016188763],
                         [0.256081040, 0.179641350, 0.06652327, 0.256081040]],
                        dtype=str(mcl.trans_matrix.dtype))
    assert np.allclose(mcl.trans_matrix, expected), print(mcl.trans_matrix)

    mcl.mcl_step()
    expected = np.array([[0.1263697, 0.1054491, 0.067736246, 0.1263697],
                         [0.7429621, 0.78150123, 0.84387976, 0.7429621],
                         [0.0042984216, 0.0076005417, 0.020647783, 0.0042984216],
                         [0.1263697, 0.1054491, 0.067736246, 0.1263697]],
                        dtype=str(mcl.trans_matrix.dtype))
    assert np.array_equal(mcl.trans_matrix, expected), print(mcl.trans_matrix)

    mcl.mcl_step()
    expected = np.array([[1.9703692e-02, 1.9275241e-02, 1.8409636e-02, 1.9703692e-02],
                         [9.6051764e-01, 9.6137083e-01, 9.6309304e-01, 9.6051764e-01],
                         [7.5001088e-05, 7.8738136e-05, 8.7761167e-05, 7.5001088e-05],
                         [1.9703692e-02, 1.9275241e-02, 1.8409636e-02, 1.9703692e-02]],
                        dtype=str(mcl.trans_matrix.dtype))
    assert np.array_equal(mcl.trans_matrix, expected), print(mcl.trans_matrix)


def test_markov_clustering_run(monkeypatch, capsys):
    data = """\
Bab\tCfu\t0.3
Bab\tOma\t0.5
Bab\tMle\t0
Cfu\tMle\t0.7
Cfu\tOma\t0.7
Oma\tMle\t0"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    mcl.run()

    expected = np.array([[0., 0., 0., 0.],
                         [1., 1., 1., 1.],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 0.]])
    assert np.array_equal(mcl.trans_matrix, expected), print(str(mcl.trans_matrix))
    mcl.clusters = [sorted(row) for row in mcl.clusters]
    assert mcl.clusters == [["Bab", "Cfu", "Mle", "Oma"]]

    data = """\
Bab\tCfu\t0.9
Bab\tOma\t0.1
Bab\tMle\t0.1
Cfu\tMle\t0.1
Cfu\tOma\t0.1
Oma\tMle\t0.9"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    mcl.run()
    expected = np.array([[1., 1., 0., 0.],
                         [0., 0., 0., 0.],
                         [0., 0., 1., 1.],
                         [0., 0., 0., 0.]])
    assert np.array_equal(mcl.trans_matrix, expected), print(str(mcl.trans_matrix))
    mcl.clusters = [sorted(row) for row in mcl.clusters]
    assert mcl.clusters == [["Mle", "Oma"], ['Bab', "Cfu"]]

    def safetyvalve_init(self, *_, **__):
        self.counter = 0
        self.global_reps = 2

    monkeypatch.setattr(br.SafetyValve, "__init__", safetyvalve_init)
    mcl = helpers.MarkovClustering(df, 2)
    capsys.readouterr()
    mcl.run()
    mcl.clusters = [sorted(row) for row in mcl.clusters]
    assert mcl.clusters == [["Mle", "Oma"], ['Bab', "Cfu"]]


def test_markov_clustering_write():
    data = """\
Bab\tCfu\t0.3
Bab\tOma\t0.5
Bab\tMle\t0
Cfu\tMle\t0.7
Cfu\tOma\t0.7
Oma\tMle\t0"""
    df = pd.read_csv(StringIO(data), sep="\t", header=None, index_col=False, names=["seq1", "seq2", "score"],
                     dtype={"seq1": "category", "seq2": "category", "score": "float32"})

    mcl = helpers.MarkovClustering(df, 2)
    mcl.run()

    tmp_file = br.TempFile()
    mcl.write(tmp_file.path)
    assert tmp_file.read() == "Bab	Cfu	Mle	Oma\n"
