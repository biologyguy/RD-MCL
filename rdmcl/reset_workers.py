#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Sep 3 2017

"""
DESCRIPTION OF PROGRAM
"""

import os
import shutil
try:
    from . import helpers
except ImportError:
    import helpers


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="reset_workers", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("wrkdb", action="store", nargs="?", help="Specify the working directory", default=os.getcwd())

    in_args = parser.parse_args()

    work_db = os.path.join(in_args.wrkdb, "work_db.sqlite")
    work_output_dir = os.path.join(in_args.wrkdb, ".worker_output")

    if os.path.isfile(work_db) and os.path.isdir(work_output_dir):
        with helpers.ExclusiveConnect(work_db) as cursor:
            root, dirs, files = next(os.walk(work_output_dir))
            for _dir in dirs:
                shutil.rmtree(os.path.join(root, _dir))
            for _file in files:
                os.remove(os.path.join(root, _file))
            cursor.execute("DELETE FROM queue")
            cursor.execute("DELETE FROM processing")
            cursor.execute("DELETE FROM complete")
            cursor.execute("DELETE FROM waiting")

            with open(os.path.join(in_args.wrkdb, "MasterClear"), "w") as ofile:
                ofile.write("240")
