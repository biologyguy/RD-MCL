#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
# from .. import rdmcl
from .. import helpers
from .. import launch_worker
import os
import sqlite3
import pandas as pd
from collections import OrderedDict
from buddysuite import buddy_resources as br
from buddysuite import AlignBuddy as Alb
# from copy import copy, deepcopy
import time
import argparse


def mock_valueerror(*args, **kwargs):
    raise ValueError(args, kwargs)


def mock_keyboardinterupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)


def test_instantiate_worker():
    temp_dir = br.TempDir()
    worker = launch_worker.Worker(temp_dir.path)
    assert worker.wrkdb_path == os.path.join(temp_dir.path, "work_db.sqlite")
    assert worker.hbdb_path == os.path.join(temp_dir.path, "heartbeat_db.sqlite")
    assert worker.output == os.path.join(temp_dir.path, ".worker_output")
    assert worker.masterclear_path == os.path.join(temp_dir.path, "MasterClear")
    assert worker.heartrate == 60
    assert type(worker.heartbeat) == launch_worker.rdmcl.HeartBeat
    assert worker.heartbeat.hbdb_path == worker.hbdb_path
    assert worker.heartbeat.pulse_rate == worker.heartrate
    assert worker.heartbeat.thread_type == "worker"
    assert worker.max_wait == 120
    assert worker.cpus == br.cpu_count() - 1
    assert worker.data_file == ""
    assert worker.start_time
    assert worker.split_time == 0
    assert worker.idle == 1
    assert worker.running == 1
    assert worker.last_heartbeat_from_master == 0
    assert worker.subjob_num == 1
    assert worker.num_subjobs == 1
    assert worker.job_id_hash is None
    worker.heartbeat.end()


def test_start_worker_no_master(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path, heartrate=1, max_wait=1)
    with pytest.raises(SystemExit):
        worker.start()
    assert worker.split_time != 0
    assert worker.last_heartbeat_from_master != 0
    worker.data_file = ".Worker_1.dat"
    out, err = capsys.readouterr()
    assert "Starting Worker_2" in out


def test_worker_check_master(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    worker.last_heartbeat_from_master = time.time() + 100
    assert worker.check_masters(20) is None

    worker.last_heartbeat_from_master = time.time() - 200
    with pytest.raises(SystemExit):
        worker.check_masters(20)
    out, err = capsys.readouterr()
    assert "Terminating Worker_None because of 2 min, 0 sec of master inactivity (spent 20% time idle)" in out


def test_worker_masterclear(hf, capsys):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)
    masterclear = temp_dir.subfile("MasterClear")
    assert os.path.isfile(worker.masterclear_path), print(worker.masterclear_path)

    with open(masterclear, "w") as ofile:
        ofile.write("12345")
    assert worker.masterclear() == 12345
    launch_worker.sys.stdout.flush()
    out, err = capsys.readouterr()
    assert out == "\r\rMasterClear signal 12345\n"
    assert not os.path.isfile(masterclear)

    with open(masterclear, "w") as ofile:
        ofile.write("Foo")
    assert worker.masterclear() == worker.max_wait
    out, err = capsys.readouterr()
    assert out == "\r\rinvalid literal for int() with base 10: 'Foo'\n"


def test_worker_clean_dead_threads(hf):
    temp_dir = br.TempDir()
    temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    worker = launch_worker.Worker(temp_dir.path)

    hb_con = sqlite3.connect(os.path.join(temp_dir.path, "heartbeat_db.sqlite"))
    hb_cursor = hb_con.cursor()
    hb_cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('master', %s)" % round(time.time()))
    master_id = hb_cursor.lastrowid
    hb_cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('worker', %s)" % round(time.time()))
    worker_id = hb_cursor.lastrowid
    hb_con.commit()

    work_con = sqlite3.connect(os.path.join(temp_dir.path, "work_db.sqlite"))
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO "
                        "waiting (hash, master_id) "
                        "VALUES ('foo', %s)" % master_id)
    work_cursor.execute("INSERT INTO "
                        "processing (hash, worker_id, master_id) "
                        "VALUES ('foo', ?, ?)", (worker_id, master_id,))
    # Orphan
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('baz', 101, 100)")
    work_con.commit()

    worker.clean_dead_threads(100)
    assert hb_cursor.execute("SELECT * FROM heartbeat WHERE thread_id=%s" % master_id).fetchone()
    assert work_cursor.execute("SELECT * FROM waiting WHERE master_id=%s" % master_id).fetchone()
    assert work_cursor.execute("SELECT * FROM processing WHERE worker_id=%s" % worker_id).fetchone()

    # No one waiting on complete
    work_cursor.execute("INSERT INTO "
                        "complete (hash, worker_id, master_id) "
                        "VALUES ('bar', ?, ?)", (worker_id, master_id,))
    work_con.commit()

    worker.clean_dead_threads(-100)
    assert not hb_cursor.execute("SELECT * FROM heartbeat WHERE thread_id=%s" % master_id).fetchone()
    assert not work_cursor.execute("SELECT * FROM waiting WHERE master_id=%s" % master_id).fetchone()
    assert not work_cursor.execute("SELECT * FROM processing WHERE worker_id=%s" % worker_id).fetchone()


def test_score_sequences(hf):
    outfile = br.TempFile()
    alb_obj = hf.get_data("cteno_panxs_aln")
    # Grab subset of alignment for manual calculation
    # Bfo-PanxαA   SQMWSQ--DDA
    # Bfr-PanxαD   V--RQIVVGGP
    alb_obj = Alb.extract_regions(alb_obj, "105:115")

    ss2_dfs = hf.get_data("ss2_dfs")
    ss2_dfs = {"Bfo-PanxαA": ss2_dfs["Bfo-PanxαA"], "Bfr-PanxαD": ss2_dfs["Bfr-PanxαD"]}
    # Update the ss2 dfs according to the alignment subsequences extracted
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
47     1  S  H      0.034       0.966       0.003
48     2  Q  H      0.071       0.926       0.004
49     3  M  H      0.371       0.649       0.003
50     4  W  C      0.802       0.211       0.004
51     5  S  C      0.852       0.151       0.010
52     6  Q  C      0.765       0.253       0.009
53     9  D  C      0.733       0.283       0.011
54    10  D  C      0.890       0.126       0.014
55    11  A  C      0.914       0.085       0.028
    '''
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
44     1  V  E      0.136       0.196       0.556
45     4  R  E      0.178       0.357       0.553
46     5  Q  H      0.157       0.530       0.525
47     6  I  H      0.114       0.771       0.319
48     7  V  H      0.217       0.675       0.212
49     8  V  C      0.461       0.443       0.166
50     9  G  C      0.837       0.040       0.077
51    10  G  C      0.940       0.008       0.063
52    11  P  C      0.606       0.015       0.402
    '''

    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].iloc[47:56]
    for indx, new in [(47, 1), (48, 2), (49, 3), (50, 4), (51, 5), (52, 6), (53, 9), (54, 10), (55, 11)]:
        ss2_dfs["Bfo-PanxαA"].set_value(indx, "indx", new)

    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].iloc[44:53]
    for indx, new in [(44, 1), (45, 4), (46, 5), (47, 6), (48, 7), (49, 8), (50, 9), (51, 10), (52, 11)]:
        ss2_dfs["Bfr-PanxαD"].set_value(indx, "indx", new)

    gap_open = -5
    gap_extend = 0

    # For score, subsmat = -0.363
    launch_worker.score_sequences([("Bfo-PanxαA", "Bfr-PanxαD", ss2_dfs["Bfo-PanxαA"], ss2_dfs["Bfr-PanxαD"])],
                                  [alb_obj, gap_open, gap_extend, outfile.path])

    assert outfile.read() == "\nBfo-PanxαA,Bfr-PanxαD,-0.3627272727272728,0.4183636363636363"


# #########  User Interface  ########## #
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

in_args = parser.parse_args([])


def test_argparse_init(monkeypatch):
    out_dir = br.TempDir()
    argv = ['launch_worker.py', '--workdb', out_dir.path]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    temp_in_args = launch_worker.argparse_init()
    assert temp_in_args.workdb == out_dir.path
    assert temp_in_args.heart_rate == 60
    assert temp_in_args.max_wait == 120
    assert not temp_in_args.log
    assert not temp_in_args.quiet


def test_main(monkeypatch, capsys):
    out_dir = br.TempDir()
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_1 because of 1 sec of master inactivity' in out

    workdb_con = sqlite3.connect(os.path.join(out_dir.path, "work_db.sqlite"))
    workdb_cursor = workdb_con.cursor()

    tables = {'queue': [(0, 'hash', 'TEXT', 0, None, 1),
                        (1, 'psi_pred_dir', 'TEXT', 0, None, 0),
                        (2, 'master_id', 'INTEGER', 0, None, 0),
                        (3, 'align_m', 'TEXT', 0, None, 0),
                        (4, 'align_p', 'TEXT', 0, None, 0),
                        (5, 'trimal', 'TEXT', 0, None, 0),
                        (6, 'gap_open', 'FLOAT', 0, None, 0),
                        (7, 'gap_extend', 'FLOAT', 0, None, 0)],
              'processing': [(0, 'hash', 'TEXT', 0, None, 1),
                             (1, 'worker_id', 'INTEGER', 0, None, 0),
                             (2, 'master_id', 'INTEGER', 0, None, 0)],
              'complete': [(0, 'hash', 'TEXT', 0, None, 1),
                           (1, 'worker_id', 'INTEGER', 0, None, 0),
                           (2, 'master_id', 'INTEGER', 0, None, 0)],
              'waiting': [(0, 'hash', 'TEXT', 0, None, 0),
                          (1, 'master_id', 'INTEGER', 0, None, 0)]}

    workdb_tables = workdb_cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    for table in tables:
        assert (table,) in workdb_tables

    for table, fields in tables.items():
        fields_query = workdb_cursor.execute("PRAGMA table_info(%s)" % table).fetchall()
        assert fields_query == fields

    hb_db_con = sqlite3.connect(os.path.join(out_dir.path, "heartbeat_db.sqlite"))
    hb_db_cursor = hb_db_con.cursor()
    tables = {'heartbeat': [(0, 'thread_id', 'INTEGER', 0, None, 1),
                            (1, 'thread_type', 'TEXT', 0, None, 0),
                            (2, 'pulse', 'INTEGER', 0, None, 0)]}

    hb_db_tables = hb_db_cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    for table in tables:
        assert (table,) in hb_db_tables

    for table, fields in tables.items():
        fields_query = hb_db_cursor.execute("PRAGMA table_info(%s)" % table).fetchall()
        assert fields_query == fields

    # Test quiet
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1", "--quiet"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    capsys.readouterr()
    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert not (out + err)

    out, err = capsys.readouterr()
    assert not (out + err)

    # Test logging mode (line breaks are inserted between start and termination messages)
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "2", "--log"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    capsys.readouterr()
    with pytest.raises(SystemExit):
        launch_worker.main()

    out, err = capsys.readouterr()
    assert "Starting Worker_3\n" in out and "Terminating Worker_3 because of 2 sec of master inactivity" in out

    # Test termination types
    monkeypatch.setattr(launch_worker.helpers, "dummy_func", mock_valueerror)
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_7 because of too many Worker crashes' in out

    monkeypatch.setattr(launch_worker.helpers, "dummy_func", mock_keyboardinterupt)

    with pytest.raises(SystemExit):
        launch_worker.main()
    out, err = capsys.readouterr()
    assert 'Terminating Worker_8 because of KeyboardInterrupt' in out
