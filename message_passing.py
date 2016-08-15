#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 11 2016
# Good SQLite tutorial: http://sebastianraschka.com/Articles/2014_sqlite_in_python_tutorial.html
import sqlite3
import os
import re
from MyFuncs import pretty_time, DynamicPrint, usable_cpu_count
import multiprocessing
from multiprocessing import cpu_count, Process, SimpleQueue, Pipe
import sys
import random
from time import time, sleep
import string


def broker_func(queue):
    sqlite_file = "sqlite_db.sqlite"
    if os.path.exists(sqlite_file):
        os.remove(sqlite_file)

    connection = sqlite3.connect(sqlite_file)
    cursor = connection.cursor()

    cursor.execute("CREATE TABLE data_table (hash TEXT PRIMARY KEY, species TEXT, gene TEXT, alignment TEXT, "
                   "graph TEXT, cluster_score REAL)")
    while True:
        if not queue.empty():
            output = queue.get()
            if output[0] == 'put':
                hash_id, field, data = output[1:4]
                cursor.execute("INSERT INTO data_table (hash, {0}) VALUES ('{1}', {2})".format(field, hash_id, data))
                connection.commit()
            elif output[0] == 'fetch':
                hash_id, field, pipe = output[1:4]
                cursor.execute("SELECT ({0}) FROM data_table WHERE hash='{1}'".format(field, hash_id))
                response = cursor.fetchone()[0]
                pipe.send(response)
            else:
                raise RuntimeError("Invalid instruction for broker. Must be either put or fetch, not %s." % output[0])

def database_populator(taxid, queue):
    queue = queue[0]
    queue.put(('put', taxid, 'cluster_score', random.random()))
    recvpipe, sendpipe = Pipe(False)
    queue.put(('fetch', taxid, 'cluster_score', sendpipe))
    response = recvpipe.recv()
    print(response)

def run_multicore_function(iterable, function, func_args=False, max_processes=0, quiet=False, out_type=sys.stdout):
    d_print = DynamicPrint(out_type)

    if max_processes == 0:
        max_processes = usable_cpu_count()

    else:
        cpus = cpu_count()
        if max_processes > cpus:
            max_processes = cpus
        elif max_processes < 1:
            max_processes = 1

    max_processes = max_processes if max_processes < len(iterable) else len(iterable)

    running_processes = 0
    child_list = []
    queue = SimpleQueue()  # Queue of tuples (hash_id, field, data)
    broker = Process(target=broker_func, args=[queue])
    broker.start()
    start_time = round(time())
    elapsed = 0
    counter = 0
    if not quiet:
        d_print.write("Running function %s() on %s cores\n" % (function.__name__, max_processes))
    # fire up the multi-core!!
    if not quiet:
        d_print.write("\tJob 0 of %s" % len(iterable))

    for next_iter in iterable:
        if type(iterable) is dict:
            next_iter = iterable[next_iter]
        while 1:     # Only fork a new process when there is a free processor.
            if running_processes < max_processes:
                # Start new process
                if not quiet:
                    d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                if func_args:
                    if not isinstance(func_args, list):
                        exit("Error in run_multicore_function(): The arguments passed into the multi-thread "
                             "function must be provided as a list")
                    func_args.insert(0, queue)
                    p = Process(target=function, args=(next_iter, func_args))

                else:
                    func_args = [queue]
                    p = Process(target=function, args=(next_iter, func_args))
                p.start()
                child_list.append(p)
                running_processes += 1
                counter += 1
                break
            else:
                # processor wait loop
                while 1:
                    for i in range(len(child_list)):
                        if child_list[i].is_alive():
                            continue
                        else:
                            child_list.pop(i)
                            running_processes -= 1
                            break

                    if not quiet:
                        if (start_time + elapsed) < round(time()):
                            elapsed = round(time()) - start_time
                            d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

                    if running_processes < max_processes:
                        break

    # wait for remaining processes to complete --> this is the same code as the processor wait loop above
    if not quiet:
        d_print.write("\tJob %s of %s (%s)" % (counter, len(iterable), pretty_time(elapsed)))

    while len(child_list) > 0:
        for i in range(len(child_list)):
            if child_list[i].is_alive():
                continue
            else:
                child_list.pop(i)
                running_processes -= 1
                break  # need to break out of the for-loop, because the child_list index is changed by pop

        if not quiet:
            if (start_time + elapsed) < round(time()):
                elapsed = round(time()) - start_time
                d_print.write("\t%s total jobs (%s, %s jobs remaining)" % (len(iterable), pretty_time(elapsed),
                                                                           len(child_list)))
    if queue.empty():
        broker.terminate()

    if not quiet:
        d_print.write("\tDONE: %s jobs in %s\n" % (len(iterable), pretty_time(elapsed)))
    # func_args = []  # This may be necessary because of weirdness in assignment of incoming arguments
    return

hashes = []

for x in range(10):
    while True:
        new_hash = random.choice(string.ascii_uppercase)
        new_hash += "".join([random.choice(string.ascii_lowercase) for _ in range(5)])
        print(new_hash)
        if new_hash in hashes:
            continue
        else:
            hashes.append(new_hash)
            break

run_multicore_function(hashes, database_populator)
