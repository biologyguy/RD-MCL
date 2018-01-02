#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Sep 3 2017

"""
DESCRIPTION OF PROGRAM
"""

import os
import shutil
from buddysuite import buddy_resources as br
try:
    from . import helpers
except ImportError:
    import helpers


def main():
    import argparse

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="reset_workers", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                     description='''\
\033[1mReset Workers\033[m
  Someone broke something... Time to wipe the slate!

  Delete all queued jobs. Don't worry, RD-MCL will queue them back up
  if it still really wants them to run.

\033[1mUsage\033[m:
  reset_workers <args>
''')

    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("-wdb", "--workdb", action="store", default=os.getcwd(), metavar="",
                              help="Specify the worker directory")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version="Reset workers version %s\n\n%s" % (helpers.VERSION,
                                                                                                       helpers.NOTICE))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()

    work_db = os.path.join(in_args.workdb, "work_db.sqlite")
    work_output_dir = os.path.join(in_args.workdb, ".worker_output")

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


if __name__ == '__main__':
    main()
