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


def test_argparse_init(monkeypatch, hf):
    out_dir = br.TempDir()
    argv = ['launch_worker.py', '--workdb', out_dir.path]
    monkeypatch.setattr(launch_worker.sys, "argv", argv)
    temp_in_args = launch_worker.argparse_init()
    assert temp_in_args.workdb == out_dir.path
    assert temp_in_args.heart_rate == 60
    assert temp_in_args.max_wait == 120
    assert not temp_in_args.log
    assert not temp_in_args.quiet
