#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Jan 12 2018

"""
Check how well a cluster relates to other clusters in an RD-MCL run

Only allow groups to be placed in larger groups
"""

try:
    from .compare_homolog_groups import prepare_clusters
    from . import helpers
except ImportError:
    from compare_homolog_groups import prepare_clusters
    import helpers

from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br
import re
import sys
from collections import OrderedDict
import os
import argparse
import pandas as pd

join = os.path.join

VERSION = helpers.VERSION
VERSION.name = "orthogroup_polisher"


def argparse_init():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="orthogroup_polisher", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mOrthogroup Polisher\033[m
  It's not really manual curation if a computer does it for you.

  Test how well RD-MCL-generated clusters fit into larger 
  clusters from that same run.

\033[1mUsage\033[m:
  orthogroup_polisher rdmcl_dir group_name [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("rdmcl_dir", action="store", help="Path to RD-MCL output directory")
    positional.add_argument("group_name", action="store", help="Name of group to test")

    # Optional commands
    #parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    #parser_flags.add_argument("--write", "-w", action="store", metavar="", help="Specify directory to write file(s)")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()
    if not os.path.isdir(in_args.rdmcl_dir):
        sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\t")
        sys.exit()

    for file_path in [join(in_args.rdmcl_dir, f) for f in ["final_clusters.txt", join("hmm", "rsquares_matrix.csv")]]:
        if not os.path.isfile(file_path):
            sys.stderr.write("Error: The provided RD-MCL output directory does not "
                             "contain the necessary file %s.\t" % file_path)
            sys.exit()

    r_squares = pd.read_csv(join(in_args.rdmcl_dir, "hmm", "rsquares_matrix.csv"))
    print(r_squares)


if __name__ == '__main__':
    main()