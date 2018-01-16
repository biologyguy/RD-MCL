#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Jan 12 2018

"""
Check how well a cluster relates to other clusters in an RD-MCL run

Only allow groups to be placed in larger groups
"""

try:
    from .compare_homolog_groups import prepare_clusters, Cluster
    from . import helpers
    #from . import rdmcl
except ImportError:
    from compare_homolog_groups import prepare_clusters, Cluster
    import helpers
    #import rdmcl

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
RED = "\033[91m"
GREEN = "\033[92m"
DEF_FONT = "\033[39m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
END = '\033[0m'


def create_truncnorm(mu, sigma, lower=0, upper=1):
    sigma = sigma if sigma > 0.001 else 0.001  # This prevents unrealistically small differences and DivBy0 errors
    dist = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    return dist


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t._ppf((1 + confidence) / 2, n - 1)
    return m, m - h, m + h


class Check(object):
    def __init__(self, rdmcl_dir):
        self.group_name = None
        self.output = None
        self.rdmcl_dir = rdmcl_dir
        self.clusters = prepare_clusters(join(self.rdmcl_dir, "final_clusters.txt"), hierarchy=True)
        self.master_clust = Cluster([seq for group, ids in self.clusters.items() for seq in ids])
        self.r_squares = pd.read_csv(join(self.rdmcl_dir, "hmm", "rsquares_matrix.csv"))
        self.within_group_rsquares = self._prepare_within_group_df()
        self.within_group_dist = create_truncnorm(helpers.mean(self.within_group_rsquares.r_square),
                                                  helpers.std(self.within_group_rsquares.r_square))
        self.between_group_rsquares = self._prepare_between_group_df()
        self.between_group_dist = create_truncnorm(helpers.mean(self.between_group_rsquares.r_square),
                                                   helpers.std(self.between_group_rsquares.r_square))

    def check(self, group_name):
        self.group_name = group_name
        query = self.clusters[group_name]
        query_score = Cluster(query, parent=self.master_clust).score()

        self.output = []
        for g, seqs in self.clusters.items():
            if g == group_name or len(seqs) < len(query):
                continue
            compare = self.r_squares.loc[((self.r_squares["rec_id1"].isin(query)) &
                                          (self.r_squares["rec_id2"].isin(seqs))) |
                                         ((self.r_squares["rec_id1"].isin(seqs)) &
                                          (self.r_squares["rec_id2"].isin(query))) &
                                         (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()

            ave, std = helpers.mean(compare.r_square), helpers.std(compare.r_square)
            upper2 = ave + (std * 2)
            upper2 = 1 if upper2 > 1 else upper2
            lower2 = ave - (std * 2)
            lower2 = 0 if lower2 < 0 else lower2
            orig_clust = Cluster(self.clusters[g], parent=self.master_clust)
            new_clust = Cluster(self.clusters[g] + query, parent=self.master_clust)
            self.output.append([g,
                                round(self.within_group_dist.cdf(upper2) - self.within_group_dist.cdf(lower2), 4),
                                round(self.between_group_dist.cdf(upper2) - self.between_group_dist.cdf(lower2), 4),
                                round(orig_clust.score() + query_score, 3),
                                round(new_clust.score(), 3)])
        self.output = sorted(self.output, key=lambda x: (x[1], -x[2]), reverse=True)

    def _prepare_within_group_df(self):
        if not os.path.isfile(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv")):
            sys.stderr.write("Preparing hmm/within_group_rsquares.csv...\n")
            within_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
            for g, seqs in self.clusters.items():
                if len(seqs) < 2:
                    continue
                clust_rsquares = self.r_squares.loc[(self.r_squares["rec_id1"].isin(seqs)) &
                                                    (self.r_squares["rec_id2"].isin(seqs)) &
                                                    (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()
                within_group_rsquares = within_group_rsquares.append(clust_rsquares, ignore_index=True)
            within_group_rsquares.to_csv(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv"))
        else:
            within_group_rsquares = pd.read_csv(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv"))
        return within_group_rsquares

    def _prepare_between_group_df(self):
        if not os.path.isfile(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv")):
            sys.stderr.write("Preparing hmm/between_group_rsquares.csv...\n")
            between_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
            i = 0
            for g1, seqs1 in self.clusters.items():
                i += 1
                if len(seqs1) < 2:
                    continue
                for g2, seqs2 in list(self.clusters.items())[i:]:
                    if len(seqs2) < 2:
                        continue
                    clust_rsquares = self.r_squares.loc[((self.r_squares["rec_id1"].isin(seqs1)) &
                                                         (self.r_squares["rec_id2"].isin(seqs2))) |
                                                        ((self.r_squares["rec_id1"].isin(seqs2)) &
                                                         (self.r_squares["rec_id2"].isin(seqs1))) &
                                                        (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()
                    between_group_rsquares = between_group_rsquares.append(clust_rsquares, ignore_index=True)
            between_group_rsquares.to_csv(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv"))
        else:
            between_group_rsquares = pd.read_csv(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv"))
        return between_group_rsquares

    def __str__(self):
        out_str = "%sTesting %s%s\n" % (BOLD, self.group_name, END)
        longest_group_name = len(sorted([g[0] for g in self.output], key=lambda x: len(x), reverse=True)[0]) + 2

        out_str += "{2}{0: <{1}}R²Within  R²Btw   OrigScore  NewScore{3}\n".format("Groups", longest_group_name,
                                                                                   UNDERLINE, END)
        for line in self.output:
            test1 = GREEN if line[1] > line[2] and line[1] >= 0.05 else RED
            test2 = GREEN if line[3] < line[4] else RED
            out_str += "{0: <{5}}{6}{1: <10}{2: <8}{7}{3: <11}{4}{8}\n".format(*line, longest_group_name, test1, test2,
                                                                               DEF_FONT)
        return out_str


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
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("--merge", "-m", action="store", metavar="", help="Name of group to add to")

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
    check.check(in_args.group_name)
    print(check)


if __name__ == '__main__':
    main()
