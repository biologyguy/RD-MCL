#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .. import reset_workers
from buddysuite import buddy_resources as br
import sys
from os.path import join
import sqlite3
from .. import helpers as hlp


def test_main(capsys, monkeypatch):
    test_dir = br.TempDir()
    workdb = test_dir.subdir("work_db")
    outputdir = test_dir.subdir("work_db/.worker_output")
    output_subdir = test_dir.subdir("work_db/.worker_output/subdirectory")
    sqlite_filepath = join(workdb, "work_db.sqlite")
    filepath = join(outputdir, "random_file")
    file = open(filepath, 'w').close()

    # Create sqlite database with correct structure and one row of data
    tables = {'queue': [('fbfbe456jkhgfeg', 'psi_pred_dir1', 'alignm1', 'alignp1', 'trimal1', 'gapopen1', 'gapextend1')],
              'processing': [('fbfbe456jkhgfeg', 'worker_id1')],
              'complete': [('fbfbe456jkhgfeg',)],
              'waiting': [('fbfbe456jkhgfeg', 'master_id1')]}

    sqlite_db = sqlite3.connect(sqlite_filepath)
    c = sqlite_db.cursor()

    c.execute('CREATE TABLE queue(hash TEXT PRIMARY KEY, psi_pred_dir TEXT, align_m TEXT, align_p TEXT, trimal TEXT, '
              'gap_open TEXT, gap_extend TEXT)')
    for value in tables['queue']:
        c.execute('INSERT INTO queue(hash, psi_pred_dir, align_m, align_p, trimal, gap_open, gap_extend) '
                  'VALUES(?,?,?,?,?,?,?)', value)

    c.execute('CREATE TABLE processing(hash TEXT PRIMARY KEY, worker_id TEXT)')
    for value in tables['processing']:
        c.execute('INSERT INTO processing(hash, worker_id) '
                  'VALUES(?,?)', value)

    c.execute('CREATE TABLE complete(hash TEXT PRIMARY KEY)')
    c.execute('INSERT INTO complete(hash) VALUES ("fbfbe456jkhgfeg")')

    c.execute('CREATE TABLE waiting(hash TEXT PRIMARY KEY, master_id TEXT)')
    for value in tables['waiting']:
        c.execute('INSERT INTO waiting(hash, master_id) VALUES (?,?)', value)

    sqlite_db.commit()
    sqlite_db.close()

    # Make sure that data was inserted in database
    for table_name, value_list in tables.items():
        sqlite_db = sqlite3.connect(sqlite_filepath)
        cursor = sqlite_db.cursor()
        cursor.execute('SELECT * FROM {tn} WHERE hash="fbfbe456jkhgfeg"'.format(tn=table_name))
        all_rows = cursor.fetchall()
        assert all_rows == value_list

    class GetCursor(object):
        def __init__(self, db_path):
            self.dbpath = db_path

        def __enter__(self):
            self.conn = sqlite3.connect(self.dbpath)
            cursor = self.conn.cursor()
            return cursor

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.conn.commit()
            self.conn.close()

    argv = ['reset_workers.py', '--workdb', workdb]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(hlp, "ExclusiveConnect", GetCursor)
    reset_workers.main()

    # Test that data was successfully deleted from database
    for table_name, value_list in tables.items():
        sqlite_db = sqlite3.connect(sqlite_filepath)
        cursor = sqlite_db.cursor()
        cursor.execute('SELECT * FROM {tn} WHERE hash="fbfbe456jkhgfeg"'.format(tn=table_name))
        all_rows = cursor.fetchall()
        assert all_rows == []
