#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 22 2015 

"""
Generate a similarity metric between two homolog groups files

- Split files up into list of groups, and each group is itself a list of ids
- Iterate over the subject list, comparing each group against the groups in the query list
    - Count the number of matches between groups
    - If one or more sequences in subject group are found in query group, divide the number of matching ids by the total
    size of both the query and subject
    - Weight the score against the total size of the subject group
    - Sum the weighted scores for all matching query groups
    - Stop iterating through query once all subject ids have been identified
    - Final tally for each subject group is weighted against the total size of all groups in subject file
    - Final score is the sum of all tallies from subject/query search, between max value of 1 and min of 0
- The metric is not symmetric between subject and query, so for a true comparison run compare() in both directions and
take the average (not currently implemented).
"""

try:
    import rdmcl
    import helpers as hlp
except ImportError:
    from rdmcl import rdmcl
    from rdmcl import helpers as hlp

import re
from collections import OrderedDict

VERSION = hlp.VERSION
VERSION.name = "compare_homolog_groups"


class Comparison(object):
    def __init__(self, true_clusters, query_clusters):
        self.true_clusters = hlp.prepare_clusters(true_clusters)
        self.query_clusters = hlp.prepare_clusters(query_clusters)

        self.total_size = len([seq_id for cluster in self.true_clusters for seq_id in cluster])

        # These scores are filled in when _prepare_difference() is called
        self.precision = 0
        self.recall = 0
        self.accuracy = 0
        self.tn_rate = 0
        self.query_score = 0
        self.true_score = 0

        self.pretty_out = ""
        self._prepare_difference()

    def _prepare_difference(self):
        # Make sure all sequences are present in both sets of clusters
        true_seq_ids = set([seq_id for clust in self.true_clusters for seq_id in clust])
        query_seq_ids = set([seq_id for clust in self.query_clusters for seq_id in clust])
        unique = [seq_id for seq_id in true_seq_ids if seq_id not in query_seq_ids]
        unique += [seq_id for seq_id in query_seq_ids if seq_id not in true_seq_ids]
        if unique:
            raise ValueError("Not all sequences are present in both sets of clusters")

        # For each true cluster, find the query cluster with the most overlap
        true_to_query = []
        query_clusters = list(self.query_clusters)

        for t_indx, t_cluster in enumerate(self.true_clusters):
            intersections = [list(filter(lambda x: x in t_cluster, sublist)) for sublist in query_clusters]

            max_match = 0
            max_match_indx = None
            for query_indx, intersect in enumerate(intersections):
                # First, we're looking for the true cluster with the highest number of query sequences
                if len(intersect) > max_match:
                    max_match = len(intersect)
                    max_match_indx = query_indx
                # If there is a tie, then we want the true cluster with the highest proportion of query seqs
                elif len(intersect) == max_match and len(intersect) > 0:
                    orig_proportion = max_match / len(query_clusters[max_match_indx])
                    new_proportion = len(intersect) / len(query_clusters[query_indx])
                    if new_proportion > orig_proportion:
                        max_match = len(intersect)
                        max_match_indx = query_indx

            true_to_query.append([t_indx, max_match_indx])

        final_clusters = [[] for _ in range(len(self.query_clusters))]
        for q_indx, q_cluster in enumerate(self.query_clusters):
            t_indx = [t_indx for t_indx, max_match_indx in true_to_query if q_indx == max_match_indx]
            if not t_indx:
                for seq_id in q_cluster:
                    if re.search("Mle", seq_id):
                        final_clusters[q_indx].append("\033[91m\033[4m%s\033[24m\033[39m" % seq_id)
                    else:
                        final_clusters[q_indx].append("\033[91m%s\033[39m" % seq_id)
            else:
                t_indx = t_indx[0]
                for seq_id in q_cluster:
                    if seq_id not in self.true_clusters[t_indx]:
                        if re.search("Mle", seq_id):
                            final_clusters[q_indx].append("\033[91m\033[4m%s\033[24m\033[39m" % seq_id)
                        else:
                            final_clusters[q_indx].append("\033[91m%s\033[39m" % seq_id)
                    else:
                        if re.search("Mle", seq_id):
                            final_clusters[q_indx].append("\033[92m\033[4m%s\033[24m\033[39m" % seq_id)
                        else:
                            final_clusters[q_indx].append("\033[92m%s\033[39m" % seq_id)

        self._metrics(true_to_query)
        for q_cluster in final_clusters:
            self.pretty_out += "%s\n" % "\t".join(q_cluster)
        return

    def _metrics(self, true_to_query):
        # Setting comparison stats, as explained in Lechner M. et al., 2014 PlosONE. DOI: 10.1371/journal.pone.0105015
        # tp = True positive, fp = False positive, fn = False negative, tn = True negative
        sum_tp, sum_fp, sum_fn, sum_tn = 0, 0, 0, 0
        for t_indx, q_indx in true_to_query:
            tp = len(list(filter(lambda x: x in self.true_clusters[t_indx], self.query_clusters[q_indx])))
            fp = len(self.query_clusters[q_indx]) - tp
            fn = len([x for x in self.true_clusters[t_indx] if x not in self.query_clusters[q_indx]])
            tn = self.total_size - tp - fp - fn

            sum_tp += tp
            sum_fp += fp
            sum_fn += fn
            sum_tn += tn

        self.precision = (sum_tp / (sum_tp + sum_fp))
        self.recall = (sum_tp / (sum_tp + sum_fn))
        self.accuracy = ((sum_tp + sum_tn) / (sum_tp + sum_tn + sum_fp + sum_fn))
        self.tn_rate = (sum_tn / (sum_tn + sum_fp))

        query_parent = Cluster([_id for _ids in self.query_clusters for _id in _ids])
        self.query_score = sum([Cluster(next_set, parent=query_parent).score() for next_set in self.query_clusters])

        true_parent = Cluster([_id for _ids in self.true_clusters for _id in _ids])
        self.true_score = sum([Cluster(next_set, parent=true_parent).score() for next_set in self.true_clusters])
        return

    def score(self):
        output = "Precision:    %s%%\n" % (round(self.precision * 100, 2))
        output += "Recall:       %s%%\n" % (round(self.recall * 100, 2))
        output += "Accuracy:     %s%%\n" % (round(self.accuracy * 100, 2))
        output += "tn rate:      %s%%\n" % (round(self.tn_rate * 100, 2))
        output += "Query score:  %s\n" % round(self.query_score, 2)
        output += "True score:   %s\n" % round(self.true_score, 2)
        return output

    def __str__(self):
        return self.score()


class Cluster(rdmcl.Cluster):
    def __init__(self, seq_ids, parent=None, taxa_separator="-"):
        self.taxa_separator = taxa_separator
        self.taxa = OrderedDict()  # key = taxa id. value = list of genes coming fom that taxa.
        self.seq_ids = sorted(seq_ids)
        for next_seq_id in seq_ids:
            taxa = next_seq_id.split(taxa_separator)[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)
        self.parent = parent
        self.cluster_score = None
        if parent:
            self.max_genes_in_a_taxa = parent.max_genes_in_a_taxa
        else:
            self.max_genes_in_a_taxa = max([len(self.taxa[taxa]) for taxa in self.taxa])


def main():
    import argparse
    from buddysuite import buddy_resources as br

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="compare_homolog_groups", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mCompare homolog groups\033[m
  An RD-MCL output visualization tool
     
  Pass in a file containing the 'true' orthogroups and a file 
  containing inferred orthogroups. A colorful comparison will
  be created!
  
\033[1mUsage\033[m:
  compare_homolog_groups "/path/to/true" "/path/to/inferred" [-options]
''')

    parser.register('action', 'setup', rdmcl._SetupAction)

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("true_clusters", help="Input file 1", action="store")
    positional.add_argument("query_clusters", help="Input file 2", action="store")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("--score", "-s", action="store_true",
                              help="Output a table of similarity instead of sequence ids")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()

    try:
        comparison = Comparison(in_args.true_clusters, in_args.query_clusters)
        if in_args.score:
            print(comparison.score())
        else:
            print(comparison.pretty_out)

    except ValueError as err:
        if "Not all sequences are present" not in str(err):
            raise

        br._stderr(hlp.RED + "Error!" + hlp.END + "\n")
        br._stderr("There are differences in the sequences present in your query and true clusters files.\n")


if __name__ == '__main__':
    main()
