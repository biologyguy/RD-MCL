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

import MyFuncs
import os
import re
import sys


class Clusters(object):
    def __init__(self, path, group_split="\n", taxa_split=" "):
        with open(path, "r") as ifile:
            self.input = ifile.read()

        self.clusters = self.input.strip().split(group_split)
        self.clusters = [[y for y in x.strip().split(taxa_split)] for x in self.clusters]
        self.clusters.reverse()
        self.size = 0.
        for group in self.clusters:
            self.size += len(group)
        self.printer = MyFuncs.DynamicPrint()

    def compare(self, query_clusters, output_file):
        score = 0.
        counter = 1
        output = ""
        for subj in self.clusters:
            printer.write("Cluster %s of %s" % (counter, len(self.clusters)))
            output += "Subj: %s\n" % subj
            counter += 1
            tally = 0.
            len_subj = len(subj)

            # This is inefficient, bc it iterates from the top of query list every time... Need a way to better manage
            # search if this ever starts to be used regularly. Possibly by converting Clusters.clusters to to a set()
            for query in query_clusters.clusters:
                matches = self.num_matches(subj, query)
                if not matches:
                    continue
                else:
                    weighted_match = (matches * 2.) / (len(subj) + len(query))
                    weighted_match *= matches / len(subj)
                    output += "Query: %s\n%s\n" % (query, weighted_match)
                    tally += weighted_match
                    len_subj -= matches
                    if len_subj == 0:
                        break

            tally *= (len(subj) / self.size)
            score += tally
            output += "Score: %s\n\n" % tally

        print("")

        if output_file:
            with open(output_file, "w") as _ofile:
                _ofile.write(output)

        return score

    @staticmethod
    def num_matches(_subj, _query):
        count = 0.
        for next_item in _subj:
            if next_item in _query:
                count += 1
        return count if count > 0 else None


def write_difference(subject, query):
    with open(subject, "r") as ifile:
        subject = ifile.readlines()
    if subject[-1] == "\n":
        del subject[-1]

    for indx, line in enumerate(subject):
        line = re.sub("group_.*?\t", "", line)
        line = re.sub("-*[0-9]+\.[0-9]*\t", "", line)
        line = line.split()
        subject[indx] = line

    with open(query, "r") as ifile:
        query = ifile.readlines()
    if query[-1] == "\n":
        del query[-1]

    for indx, line in enumerate(query):
        line = re.sub("group_.*?\t", "", line)
        line = re.sub("-*[0-9]+\.[0-9]*\t", "", line)
        line = line.split()
        query[indx] = line

    # For each query cluster, find the subject cluster with the most overlap
    final_clusters = [[] for _ in range(len(query))]
    for query_indx, cluster in enumerate(query):
        intersections = [list(filter(lambda x: x in cluster, sublist)) for sublist in subject]
        max_match = 0
        max_match_indx = None
        for sub_indx, intersect in enumerate(intersections):
            if len(intersect) > max_match:
                max_match = len(intersect)
                max_match_indx = sub_indx

        for seq_id in cluster:
            if max_match_indx is None or seq_id not in subject[max_match_indx]:
                final_clusters[query_indx].append("\033[91m%s\033[39m" % seq_id)
            else:
                final_clusters[query_indx].append("\033[92m%s\033[39m" % seq_id)
        if max_match_indx is not None:
            subject[max_match_indx] = []

    for cluster in final_clusters:
        print("\t".join(cluster))
    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="compare_homolog_groups", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("subject", help="Input file 1", action="store")
    parser.add_argument("query", help="Input file 2", action="store")
    parser.add_argument("--group_split", "-gs", action="store", default="\n",
                        help="specify the delimiting string between groups")
    parser.add_argument("--taxa_split", "-ts", action="store", default=" ",
                        help="specify the delimiting string between taxa")
    parser.add_argument("--output_file", "-o", help="Specify a location to send an output file.", action="store")

    in_args = parser.parse_args()

    timer = MyFuncs.Timer()
    printer = MyFuncs.DynamicPrint()

    ofile = None if not in_args.output_file else os.path.abspath(in_args.output_file)

    #groups1 = Clusters(in_args.subject, in_args.group_split, in_args.taxa_split)
    #groups2 = Clusters(in_args.query, in_args.group_split, in_args.taxa_split)

    #print("Score: %s\n%s" % (groups1.compare(groups2, ofile), timer.total_elapsed()))

    write_difference(in_args.subject, in_args.query)