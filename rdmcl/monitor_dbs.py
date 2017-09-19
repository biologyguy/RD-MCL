#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Jul 25 2017 

"""
DESCRIPTION OF PROGRAM
"""

import os
from multiprocessing import Process
from time import time, sleep
from buddysuite.buddy_resources import DynamicPrint, TempFile

try:
    import helpers
except ImportError:
    from . import helpers


class Monitor(object):
    def __init__(self, hbdb_path, wdb_path):
        self.hbdb_path = hbdb_path
        self.wdb_path = wdb_path

    def _run(self, check_file_path):
        printer = DynamicPrint()
        output = [("#Master", 9), ("AveMhb", 9), ("#Worker", 9), ("AveWhb", 9), ("#queue", 9),
                  ("#subq", 8), ("#proc", 8), ("#subp", 8), ("#comp", 8), ("#HashWait", 10),
                  ("#IdWait", 9), ("ConnectTime", 9)]

        output = [str(x[0]).ljust(x[1]) for x in output]
        printer.write("".join(output))
        printer.new_line(1)
        while True:
            with open("%s" % check_file_path, "r") as ifile:
                ifile_content = ifile.read()
            if ifile_content != "Running":
                break

            split_time = time()
            try:
                with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                    heartbeat = cursor.execute("SELECT * FROM heartbeat").fetchall()
                with helpers.ExclusiveConnect(self.wdb_path) as cursor:
                    queue = cursor.execute("SELECT hash FROM queue").fetchall()
                    processing = cursor.execute("SELECT hash FROM processing").fetchall()
                    complete = cursor.execute("SELECT master_id, worker_id FROM complete").fetchall()
                    waiting = cursor.execute("SELECT * FROM waiting").fetchall()
                connect_time = round(time() - split_time, 3)
                master_hb = [hb[2] for hb in heartbeat if hb[1] == "master"]
                len_mas = len(master_hb)
                master_hb = 0 if not master_hb else round(time() - (sum(master_hb) / len(master_hb)), 1)
                worker_hb = [hb[2] for hb in heartbeat if hb[1] == "worker"]
                len_wor = len(worker_hb)
                worker_hb = 0 if not worker_hb else round(time() - (sum(worker_hb) / len(worker_hb)), 1)

                hash_wait = {}
                master_wait = {}
                for seq_hash, master_id in waiting:
                    hash_wait.setdefault(seq_hash, 0)
                    hash_wait[seq_hash] += 1
                    master_wait.setdefault(master_id, 0)
                    master_wait[master_id] += 1

                subqueue_len = len([None for x in queue if "_" in x[0]])
                queue_len = len(queue) - subqueue_len

                subproc_len = len([None for x in processing if "_" in x[0]])
                proc_len = len(processing) - subproc_len

                output = [(len_mas, 9), (master_hb, 9), (len_wor, 9), (worker_hb, 9), (queue_len, 9),
                          (subqueue_len, 8), (proc_len, 8), (subproc_len, 8), (len(complete), 8),
                          (len(hash_wait), 10), (len(master_wait), 9), (connect_time, 9)]
                output = [str(x[0]).ljust(x[1]) for x in output]

                printer.write("".join(output))
            except helpers.sqlite3.OperationalError:
                printer.write("Worker databases not detected")
            sleep(0.5)
        return

    def start(self):
        tmp_file = TempFile()
        tmp_file.write("Running")
        p = Process(target=self._run, args=(tmp_file.path,))
        p.daemon = 1
        p.start()
        input("Press return to terminate.\n")
        tmp_file.clear()
        p.join()
        return


def main():
    import argparse

    parser = argparse.ArgumentParser(prog="monitor_dbs", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(),
                        help="Specify the directory where sqlite databases will be fed by RD-MCL", )

    in_args = parser.parse_args()

    workdb = os.path.join(in_args.workdb, "work_db.sqlite")
    heartbeatdb = os.path.join(in_args.workdb, "heartbeat_db.sqlite")

    monitor = Monitor(heartbeatdb, workdb)
    monitor.start()


if __name__ == '__main__':
    main()
