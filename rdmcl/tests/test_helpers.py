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


def test_sqlitebroker_init(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    assert broker.db_file == "%s%sdb.sqlite" % (tmpdir.path, hf.sep)
    assert type(broker.connection) == sqlite3.Connection
    assert type(broker.cursor) == sqlite3.Cursor
    assert type(broker.broker_queue) == SimpleQueue
    assert broker.broker is None
