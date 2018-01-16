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

#from buddysuite import SeqBuddy as Sb
#from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br
#import re
import sys
#from collections import OrderedDict
import os
import argparse
import pandas as pd
import numpy as np
from scipy import stats

join = os.path.join

VERSION = helpers.VERSION
VERSION.name = "orthogroup_polisher"


def create_truncnorm(mu, sigma, lower=0, upper=1):
    sigma = sigma if sigma > 0.001 else 0.001  # This prevents unrealistically small differences and DivBy0 errors
    dist = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    return dist


def group_compare(group1, group2, df):
    compared_rsquares_df = df.loc[((df["rec_id1"].isin(group1)) &
                                   (df["rec_id2"].isin(group2)))|
                                  ((df["rec_id1"].isin(group2)) &
                                   (df["rec_id2"].isin(group1)))&
                                  (df["rec_id1"] != df["rec_id2"])].copy()
    return compared_rsquares_df


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


def check_group(group_name, clusts, within_distr, btw_distr, all_r_squares):
    query = clusts[group_name]
    print("Testing", group_name)
    print("Group\t\tWithin\tBetween")
    output = []
    for g, seqs in clusts.items():
        if g == group_name or len(seqs) < len(query):
            continue
        compare = group_compare(query, seqs, all_r_squares)
        ave, std = helpers.mean(compare.r_square), helpers.std(compare.r_square)
        upper2 = ave + (std * 2)
        upper2 = 1 if upper2 > 1 else upper2
        lower2 = ave - (std * 2)
        lower2 = 0 if lower2 < 0 else lower2
        output.append([g,
                       round(within_distr.cdf(upper2) - within_distr.cdf(lower2), 5),
                       round(btw_distr.cdf(upper2) - btw_distr.cdf(lower2), 5)])
    output = sorted(output, key=lambda x: (x[1], -x[2]), reverse=True)
    for line in output:
        print("\t".join([str(l) for l in line]))


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
    rdmcl_dir = os.path.abspath(in_args.rdmcl_dir)
    if not os.path.isdir(rdmcl_dir):
        sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\n")
        sys.exit()

    if not os.path.isfile(join(rdmcl_dir, "final_clusters.txt")):
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'final_clusters.txt'.\n")
        sys.exit()
    clusters = prepare_clusters(join(rdmcl_dir, "final_clusters.txt"), hierarchy=True)

    if not os.path.isfile(join(rdmcl_dir, "hmm", "rsquares_matrix.csv")):
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'hmm/rsquares_matrix.csv'.\n")
        sys.exit()
    r_squares = pd.read_csv(join(rdmcl_dir, "hmm", "rsquares_matrix.csv"))

    if not os.path.isfile(join(rdmcl_dir, "hmm", "within_group_rsquares.csv")):
        sys.stderr.write("Preparing hmm/within_group_rsquares.csv...\n")
        within_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
        for g, seqs in clusters.items():
            if len(seqs) < 2:
                continue
            clust_rsquares = r_squares.loc[(r_squares["rec_id1"].isin(seqs)) &
                                           (r_squares["rec_id2"].isin(seqs)) &
                                           (r_squares["rec_id1"] != r_squares["rec_id2"])].copy()
            within_group_rsquares = within_group_rsquares.append(clust_rsquares, ignore_index=True)
        within_group_rsquares.to_csv(join(rdmcl_dir, "hmm", "within_group_rsquares.csv"))
    else:
        within_group_rsquares = pd.read_csv(join(rdmcl_dir, "hmm", "within_group_rsquares.csv"))

    if not os.path.isfile(join(rdmcl_dir, "hmm", "between_group_rsquares.csv")):
        sys.stderr.write("Preparing hmm/between_group_rsquares.csv...\n")
        between_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
        i = 0
        for g1, seqs1 in clusters.items():
            i += 1
            if len(seqs1) < 2:
                continue
            for g2, seqs2 in list(clusters.items())[i:]:
                if len(seqs2) < 2:
                    continue
                clust_rsquares = r_squares.loc[((r_squares["rec_id1"].isin(seqs1)) &
                                                (r_squares["rec_id2"].isin(seqs2))) |
                                               ((r_squares["rec_id1"].isin(seqs2)) &
                                                (r_squares["rec_id2"].isin(seqs1))) &
                                               (r_squares["rec_id1"] != r_squares["rec_id2"])].copy()
                between_group_rsquares = between_group_rsquares.append(clust_rsquares, ignore_index=True)
        between_group_rsquares.to_csv(join(rdmcl_dir, "hmm", "between_group_rsquares.csv"))
    else:
        between_group_rsquares = pd.read_csv(join(rdmcl_dir, "hmm", "between_group_rsquares.csv"))

    # Pretend that each distribution is normal, even though they are complex
    within_group_dist = create_truncnorm(helpers.mean(within_group_rsquares.r_square),
                                         helpers.std(within_group_rsquares.r_square))
    between_group_dist = create_truncnorm(helpers.mean(between_group_rsquares.r_square),
                                          helpers.std(between_group_rsquares.r_square))

    check_group(in_args.group_name, clusters, within_group_dist, between_group_dist, r_squares)


if __name__ == '__main__':
    main()