#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 11 2016
# Good SQLite tutorial: http://sebastianraschka.com/Articles/2014_sqlite_in_python_tutorial.html
import sqlite3
import os
import re

psipred_dir = "rdmcl/tmp_dir/psi_pred"

'''
Connect to sqlite database (file can be stored anywhere on disk)
Cursor executes commands on the database
'''
sqlite_file = "sqlite_db.sqlite"
if os.path.exists(sqlite_file):
    os.remove(sqlite_file)

connection = sqlite3.connect(sqlite_file)
cursor = connection.cursor()


'''
Data types supported:
INTEGER -> integer up to 8 bytes
REAL -> float up to 8 bytes
TEXT -> (typically) utf-8 string
BLOB -> binary data
NULL -> missing data/empty cell
No boolean type, use a 1 or 0 with the INTEGER type
'''

'''
Common Operations (called within cursor.execute(<command>)

Table creation/data insertion:
CREATE TABLE <table_name> (<column_name> <column_type>, ..., <optional>)
    Adds a new table with as many fields as desired.
    Note: Can add PRIMARY KEY to the end of a column declaration to specify it as the 'primary key'
    Using a primary key to search is much faster than other columns, since it is guaranteed to be UNIQUE
    CONSTRAINT <primary_key_column> PRIMARY KEY (<column>, ...)
        Can set the primary key as one or more columns (combined columns)
ALTER TABLE <table_name> <sub_command> -> Alters the table.
    ADD COLUMN '<column_name>' <column_type> <optional>
    Adds a column to the table.
        DEFAULT <value> (optional)
        Sets the default value. If not provided, value is set to NULL.
INSERT INTO <table_name> (<column_name>, ...) VALUES (<value>, ...)
    Adds new columns to the table with values
    Note: If the one of the items is being added to the primary key column and already exists in another row, an
    IntegrityError will be thrown.
INSERT OR IGNORE INTO <table_name> (<column_name>, ...) VALUES (<value>, ...)
    Same as INSERT INTO but if there is an exception it ignores the command
UPDATE <table_name> SET <column_name>=('<value>') WHERE <column_to_search>=(<query>)
    Using the PRIMARY KEY as the search column is preferred for speed
CREATE INDEX <index_name> on <table_name>(<new_column>)
    Creates a column which requires UNIQUE values, for fast searching
DROP INDEX <index_name>
    Drops an index from the table

Data retrieval (should follow up with connection.fetchall():
* = all
SELECT <column_id or *>,... FROM <table_name> WHERE <column_to_search>=<query> <optional>
    Returns the value of specified columns from rows matching the query
    LIMIT <num> (optional)
    Specifies the maximum number of rows to be matched

PRAGMA TABLE_INFO(<table_name>)
    Returns a list of tuples, where each tuple represents a row. If no primary key is specified, the first index of the
    tuple will be a row index integer.

'''
cursor.execute("CREATE TABLE psipred_table (seq_id TEXT PRIMARY KEY, species TEXT, gene TEXT, data TEXT)")
cursor.execute("CREATE TABLE alignment_table (group_id TEXT PRIMARY KEY, index1 INTEGER, index2 INTEGER, data TEXT)")

psipred_files = os.listdir(psipred_dir)
for ps_filename in psipred_files:
    if '-' in ps_filename and '.ss2' in ps_filename:
        seq_name = re.sub('\.ss2', '', ps_filename)
        species, gene = seq_name.split('-')[0], seq_name.split('-')[1]
        with open("{0}/{1}".format(psipred_dir, ps_filename), 'r') as psi_data:
            data = psi_data.read()
        cursor.execute("INSERT INTO psipred_table (seq_id, species, gene, data) VALUES ('{0}', '{1}', '{2}', '{3}')"
                       .format(seq_name, species, gene, data))

cursor.execute("SELECT (gene) FROM psipred_table WHERE seq_id='Bfo-PanxαI'")
result = cursor.fetchall()
#print(result[0][0])

cursor.execute("SELECT (seq_id) FROM psipred_table WHERE species='Bfo'")
result = cursor.fetchall()
#print(result)

cursor.execute("SELECT * FROM psipred_table WHERE species='Bfo'")
result = cursor.fetchall()

cursor.execute("SELECT (seq_id) FROM psipred_table WHERE species='ABCDEFGHIJKLMNOPQRSTUVWXYZ'")
result = cursor.fetchall()


connection.commit()  # All changes must be committed to the database before they appear
connection.close()  # Must close the database before exiting the program
