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
    from rdmcl import rdmcl
except ImportError:
    import rdmcl

import os
import re
import sys
from collections import OrderedDict


class Comparison(object):
    def __init__(self, true_clusters, query_clusters):
        self.true_clusters = self.prepare_clusters(true_clusters)
        self.query_clusters = self.prepare_clusters(query_clusters)

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

    @staticmethod
    def prepare_clusters(ifile):
        with open(ifile, "r") as ifile:
            output = ifile.readlines()
        if output[-1] == "\n":
            del output[-1]

        for indx, line in enumerate(output):
            line = re.sub("group_.*?\t", "", line)
            line = re.sub("^-*[0-9]+\.[0-9]*\t", "", line)
            line = line.split()
            output[indx] = line
        return output

    def _prepare_difference(self):
        # For each query cluster, find the true cluster with the most overlap
        final_clusters = [[] for _ in range(len(self.query_clusters))]
        query_to_true = []
        success_tally = 0
        true_clusters = list(self.true_clusters)

        for q_indx, q_cluster in enumerate(self.query_clusters):
            intersections = [list(filter(lambda x: x in q_cluster, sublist)) for sublist in true_clusters]
            max_match = 0
            max_match_indx = None
            for true_indx, intersect in enumerate(intersections):
                if len(intersect) > max_match:
                    max_match = len(intersect)
                    max_match_indx = true_indx

            query_to_true.append([q_indx, max_match_indx])
            for seq_id in q_cluster:
                if max_match_indx is None or seq_id not in self.true_clusters[max_match_indx]:
                    if re.search("Mle", seq_id):
                        final_clusters[q_indx].append("\033[91m\033[4m%s\033[24m\033[39m" % seq_id)
                    else:
                        final_clusters[q_indx].append("\033[91m%s\033[39m" % seq_id)
                else:
                    success_tally += 1
                    if re.search("Mle", seq_id):
                        final_clusters[q_indx].append("\033[92m\033[4m%s\033[24m\033[39m" % seq_id)
                    else:
                        final_clusters[q_indx].append("\033[92m%s\033[39m" % seq_id)
            if max_match_indx is not None:
                true_clusters[max_match_indx] = []

        # Setting comparison stats, as explained in Lechner M. et al., 2014 PlosONE. DOI: 10.1371/journal.pone.0105015
        # tp = True positive, fp = False positive, fn = False negative, tn = True negative
        for q_indx, t_indx in query_to_true:
            if t_indx is None:
                tp = 0
                fp = len(self.query_clusters[q_indx])
                fn = 0
                tn = self.total_size - fp
            else:
                tp = len(list(filter(lambda x: x in self.query_clusters[q_indx], self.true_clusters[t_indx])))
                fp = len(self.query_clusters[q_indx]) - tp
                fn = len([x for x in self.true_clusters[t_indx] if x not in self.query_clusters[q_indx]])
                tn = self.total_size - tp - fp - fn

            self.precision += (tp / (tp + fp)) * (len(self.query_clusters[q_indx]) / self.total_size)
            self.recall += 0 if t_indx is None else (tp / (tp + fn)) * (len(self.query_clusters[q_indx]) / self.total_size)
            self.accuracy += ((tp + tn) / self.total_size) * (len(self.query_clusters[q_indx]) / self.total_size)
            self.tn_rate += (tn / (tn + fp)) * (len(self.query_clusters[q_indx]) / self.total_size)

        query_parent = Cluster([_id for _ids in self.query_clusters for _id in _ids])
        self.query_score = sum([Cluster(next_set, parent=query_parent).score() for next_set in self.query_clusters])

        true_parent = Cluster([_id for _ids in self.true_clusters for _id in _ids])
        self.true_score = sum([Cluster(next_set, parent=true_parent).score() for next_set in self.true_clusters])

        for q_cluster in final_clusters:
            self.pretty_out += "%s\n" % "\t".join(q_cluster)
        return


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


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="compare_homolog_groups", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("true_clusters", help="Input file 1", action="store")
    parser.add_argument("query_clusters", help="Input file 2", action="store")
    parser.add_argument("--score", "-s", action="store_true")
    parser.add_argument("--group_split", "-gs", action="store", default="\n",
                        help="specify the delimiting string between groups")
    parser.add_argument("--taxa_split", "-ts", action="store", default=" ",
                        help="specify the delimiting string between taxa")
    parser.add_argument("--output_file", "-o", help="Specify a location to send an output file.", action="store")

    in_args = parser.parse_args()

    ofile = None if not in_args.output_file else os.path.abspath(in_args.output_file)

    comparison = Comparison(in_args.true_clusters, in_args.query_clusters)

    if in_args.score:
        print("Precision:    %s%%" % (round(comparison.precision, 4) * 100))
        print("Recall:       %s%%" % (round(comparison.recall, 4) * 100))
        print("Accuracy:     %s%%" % (round(comparison.accuracy, 4) * 100))
        print("tn rate:      %s%%" % (round(comparison.tn_rate, 4) * 100))
        print("Query score:  %s" % round(comparison.query_score, 2))
        print("True score:   %s\n" % round(comparison.true_score, 2))
    else:
        print(comparison.pretty_out)
