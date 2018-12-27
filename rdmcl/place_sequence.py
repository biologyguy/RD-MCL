#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Jan 18 2018

"""
Given RD-MCL clusters, where do new sequences fit in?
"""

try:
    from rdmcl import rdmcl
    from rdmcl import helpers as hlp
    from rdmcl.merge_orthogroups import Check
except ImportError:
    import rdmcl
    import helpers as hlp
    from merge_orthogroups import Check

from buddysuite import buddy_resources as br
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
import sys
import os
from os.path import join
import argparse
import numpy as np

VERSION = hlp.VERSION
VERSION.name = "place_sequences"


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
    parser_flags.add_argument("--min", "-min", action="store", metavar="", type=int, default=5,
                              help="Only test groups of min size")
    parser_flags.add_argument("--force", "-f", action="store_true",
                              help="Automatically answer 'yes' to any warning messages. Use caution!")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    print("%sPlacement of sequence(s) in RD-MCL clusters%s\n" % (hlp.BOLD, hlp.END))
    in_args = argparse_init()
    rdmcl_dir = os.path.abspath(in_args.rdmcl_dir)

    if not os.path.isdir(rdmcl_dir):
        sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\n")
        sys.exit()

    check_files = ["final_clusters.txt", join("hmm", "rsquares_matrix.csv"),
                   join("hmm", "hmm_fwd_scores.csv"), "input_seqs.fa"]
    for check_file in check_files:
        if not os.path.isfile(join(rdmcl_dir, check_file)):
            sys.stderr.write("Error: The provided RD-MCL output directory does not "
                             "contain the necessary file '%s'.\n" % check_file)
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

    upper_r2_btw_groups = check.btw_group_r2_dist.mean() + (check.btw_group_r2_dist.std() * 2)
    upper_r2_btw_groups = 1. if upper_r2_btw_groups > 1 else upper_r2_btw_groups
    lower_r2_btw_groups = check.btw_group_r2_dist.mean() - (check.btw_group_r2_dist.std() * 2)
    lower_r2_btw_groups = 0. if lower_r2_btw_groups < 0. else lower_r2_btw_groups
    print("R² between clusters: {0:0<5} - {1:0<5}".format(round(lower_r2_btw_groups, 3),
                                                          round(upper_r2_btw_groups, 3)))

    upper_r2_win_groups = check.within_group_r2_dist.mean() + (check.within_group_r2_dist.std() * 2)
    upper_r2_win_groups = 1. if upper_r2_win_groups > 1 else upper_r2_win_groups
    lower_r2_win_groups = check.within_group_r2_dist.mean() - (check.within_group_r2_dist.std() * 2)
    lower_r2_win_groups = 0. if lower_r2_win_groups < 0. else lower_r2_win_groups
    print("R² within clusters: {0:0<5} - {1:0<5}\n".format(round(lower_r2_win_groups, 3),
                                                           round(upper_r2_win_groups, 3)))

    resample = check.btw_group_fwd_dist.resample(10000)[0]
    clique95 = [np.percentile(resample, 2.5), np.percentile(resample, 97.5)]

    print("FwdScore between clusters: {0:0<5} - {1:0<5}".format(round(clique95[0], 3),
                                                                round(clique95[1], 3)))

    resample = check.within_group_fwd_dist.resample(10000)[0]
    clique95 = [np.percentile(resample, 2.5), np.percentile(resample, 97.5)]

    print("FwdScore within clusters: {0:0<5} - {1:0<5}\n".format(round(clique95[0], 3),
                                                                 round(clique95[1], 3)))

    for rec in seqbuddy.records:
        check.check_new_sequence(rec, in_args.min)

        print("%sTesting %s%s: %s" % (hlp.BOLD, rec.id, hlp.END, rec.description))
        longest_group_name = len(sorted([g[0] for g in check.output], key=lambda x: len(x), reverse=True)[0]) + 2

        out_str = "{4}{0: <{1}}Size  {2: <17}{3: <17}{5}\n".format("Group", longest_group_name, "R² 95%-CI",
                                                                   "Forward 95%-CI", hlp.UNDERLINE, hlp.END)

        longest_fwd_score = len(str(check.output[0][4]))
        for output in check.output:
            out_str += "{0: <{6}}{1: <6}{2:0<6} - {3:0<6}  {4: <{7}} - {5: <7}\n".format(*output,
                                                                                         longest_group_name,
                                                                                         longest_fwd_score)
        print(out_str, "\n")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        br._stderr("\n")
