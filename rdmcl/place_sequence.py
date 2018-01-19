#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Jan 18 2018

"""
Given RD-MCL clusters, where do new sequences fit in?
"""

try:
    from .merge_orthogroups import Check
    from . import helpers as hlp
    from . import rdmcl
except ImportError:
    from merge_orthogroups import Check
    import helpers as hlp
    import rdmcl

from buddysuite import buddy_resources as br
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
import sys
import os
from os.path import join
import pandas as pd
from subprocess import Popen, PIPE
from multiprocessing import Lock, Process
from io import StringIO
import argparse

VERSION = hlp.VERSION
VERSION.name = "merge_orthogroups"
LOCK = Lock()


def fwd_back_run(self, rec, args):
    hmm_scores_file = args[0]
    hmm_path = join(self.outdir, "hmm", "%s.hmm" % rec.id)

    fwdback_output = Popen("%s %s %s" % (rdmcl.HMM_FWD_BACK, hmm_path, self.tmp_dir.subfiles[0]),
                           shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].decode()
    fwd_scores_df = pd.read_csv(StringIO(fwdback_output), delim_whitespace=True,
                                header=None, comment="#", index_col=False)
    fwd_scores_df.columns = ["rec_id", "fwd_raw", "back_raw", "fwd_bits", "back_bits"]
    fwd_scores_df["hmm_id"] = rec.id

    hmm_fwd_scores = pd.DataFrame(columns=["hmm_id", "rec_id", "fwd_raw"])
    hmm_fwd_scores = hmm_fwd_scores.append(fwd_scores_df.loc[:, ["hmm_id", "rec_id", "fwd_raw"]],
                                           ignore_index=True)
    hmm_fwd_scores = hmm_fwd_scores.to_csv(path_or_buf=None, header=None, index=False, index_label=False)
    with LOCK:
        with open(hmm_scores_file, "a") as ofile:
            ofile.write(hmm_fwd_scores)
    return


def argparse_init():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="place_sequence", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mPlace Sequence\033[m
  Oh sure, noooowwww you want to include that sequence....

  Assess whether new sequences fit in RD-MCL derived orthogroups 

\033[1mUsage\033[m:
  place_sequence sequence_file/stdin rdmcl_dir [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("sequence", action="store", help="Query sequence(s)")
    positional.add_argument("rdmcl_dir", action="store", help="Path to RD-MCL output directory")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("--merge", "-m", action="store", metavar="", help="Name of group to add to")
    parser_flags.add_argument("--force", "-f", action="store_true",
                              help="Automatically answer 'yes' to any warning messages. Use caution!")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()
    rdmcl_dir = os.path.abspath(in_args.rdmcl_dir)

    if not os.path.isdir(rdmcl_dir):
        sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\n")
        sys.exit()

    if not os.path.isfile(join(rdmcl_dir, "final_clusters.txt")):
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'final_clusters.txt'.\n")
        sys.exit()

    if not os.path.isfile(join(rdmcl_dir, "hmm", "rsquares_matrix.csv")):
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'hmm/rsquares_matrix.csv'.\n")
        sys.exit()

    check = Check(in_args.rdmcl_dir)
    for seq_id in [seq for g_name, clust in check.clusters.items() for seq in clust]:
        if not os.path.isfile(join(rdmcl_dir, "hmm", "%s.hmm" % seq_id)):
            sys.stderr.write("Error: The provided RD-MCL output directory does not "
                             "contain the necessary HMM file for %s.\n" % seq_id)
            sys.exit()

    seqbuddy = Sb.SeqBuddy(in_args.sequence)
    hmm_dir = join(in_args.rdmcl_dir, "hmm")
    for rec in seqbuddy.records:
        if not os.path.isfile(join(hmm_dir, "%s.hmm" % rec.id)):
            align = Alb.AlignBuddy(rec.format("fasta"))
            align = Alb.generate_hmm(align, rdmcl.HMMBUILD)
            with open(join(hmm_dir, "%s.hmm" % rec.id), "w") as _ofile:
                _ofile.write(align.alignments[0].hmm)

    for rec in seqbuddy.records:
        check.check_new_sequence(rec)
        print(check.output)

    if len(seqbuddy) > 1 and in_args.merge:
        sys.stderr.write("Error: --merge flag provided by multiple query seqeunces present. Can only place one "
                         "sequence at a time.")
        sys.exit()

    if in_args.merge:
        #check.merge(in_args.merge, in_args.force)
        pass

if __name__ == '__main__':
    main()
