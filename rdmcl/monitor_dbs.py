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
        output = [("#Master", 10), ("AveMhb", 10), ("#Worker", 10), ("AveWhb", 10), ("#queue", 10),
                  ("#proc", 10), ("#comp", 10), ("#HashWait", 11),
                  ("#IdWait", 10), ("ConnectTime", 10)]

        output = [str(x[0]).ljust(x[1]) for x in output]
        printer.write("".join(output))
        printer.new_line(1)
        while True:
            with open("%s" % check_file_path, "r") as ifile:
                ifile_content = ifile.read()
            if ifile_content != "Running":
                break

            split_time = time()
            if os.path.isfile(self.hbdb_path) and os.path.isfile(self.wdb_path):
                with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                    heartbeat = cursor.execute("SELECT * FROM heartbeat").fetchall()
                with helpers.ExclusiveConnect(self.wdb_path) as cursor:
                    queue = cursor.execute("SELECT master_id FROM queue").fetchall()
                    processing = cursor.execute("SELECT master_id, worker_id FROM processing").fetchall()
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

                output = [(len_mas, 10), (master_hb, 10), (len_wor, 10), (worker_hb, 10), (len(queue), 10),
                          (len(processing), 10), (len(complete), 10), (len(hash_wait), 11),
                          (len(master_wait), 10), (connect_time, 10)]
                output = [str(x[0]).ljust(x[1]) for x in output]

                printer.write("".join(output))
            else:
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


if __name__ == '__main__':
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
