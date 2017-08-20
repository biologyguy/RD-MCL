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
# from copy import copy, deepcopy
import shutil
import argparse


def mock_valueerror(*args, **kwargs):
    raise ValueError(args, kwargs)


def mock_keyboardinterupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)

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
        assert 'Terminating Worker_1 because of 1 sec of maser inactivity' in err

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

    # Test logging mode
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "2", "--log"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    capsys.readouterr()
    with pytest.raises(SystemExit):
        launch_worker.main()

    out, err = capsys.readouterr()
    assert out == """
Starting Worker_3

Terminating Worker_3 because of deleted check file.
"""

    # Test termination types
    monkeypatch.setattr(launch_worker.rdmcl.HeartBeat, "start", mock_valueerror)
    argv = ['launch_worker.py', '--workdb', out_dir.path, "--max_wait", "1"]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)

    with pytest.raises(SystemExit):
        launch_worker.main()
        out, err = capsys.readouterr()
        assert 'Terminating Worker_1 because of too many Worker crashes' in err

    monkeypatch.setattr(launch_worker.rdmcl.HeartBeat, "start", mock_keyboardinterupt)

    with pytest.raises(SystemExit):
        launch_worker.main()
        out, err = capsys.readouterr()
        assert 'Terminating Worker_1 because of KeyboardInterrupt' in err
