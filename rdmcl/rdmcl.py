#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: rdmcl.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/RD-MCL
© license: None, this work is public domain

Description:
Identify 'orthogroup' from a collection of sequences. MCMCMC-driven MCL is the basic workhorse for clustering, and then
the groups are refined further to include orphan sequences and/or be broken up into RBH-cliques, where appropriate.
"""

# Std library
import sys
import os
import re
import shutil
import json
import logging
import sqlite3
from io import StringIO
from subprocess import Popen, PIPE, check_output, CalledProcessError
from multiprocessing import Lock, Pipe
from random import gauss, Random, randint
from math import ceil, log2
from collections import OrderedDict
from copy import deepcopy

# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo

# My packages
import mcmcmc
import helpers
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br

# Globals
try:
    script_path = os.path.abspath(__file__).split(os.sep)
    script_path = os.sep + os.path.join(*script_path[:-1])
    git_commit = check_output(['git', '--git-dir={0}{1}..{1}.git'.format(script_path, os.sep), 'rev-parse',
                               '--short', 'HEAD']).decode().strip()
    git_commit = " (git %s)" % git_commit if git_commit else ""
    VERSION = "1.alpha%s" % git_commit
except CalledProcessError:
    VERSION = "1.alpha"

NOTICE = '''\
Public Domain Notice
--------------------
This is free software; see the LICENSE for further details.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov
'''
LOCK = Lock()
CPUS = br.usable_cpu_count()
TIMER = helpers.Timer()
GAP_OPEN = -5
GAP_EXTEND = 0
BLOSUM62 = helpers.make_full_mat(SeqMat(MatrixInfo.blosum62))
MCMCMC_CHAINS = 3

ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
               "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
for aa in ambiguous_X:
    pair = sorted((aa, "X"))
    pair = tuple(pair)
    BLOSUM62[pair] = ambiguous_X[aa]


def push(hash_id, field, data):
    """This is deprecated and needs to be removed"""
    broker_queue.put(('push', hash_id, field, data))


def scored_clusters():
    """This is deprecated and needs to be removed"""
    recvpipe, sendpipe = Pipe(False)
    broker_queue.put(('scored', sendpipe))
    response = json.loads(recvpipe.recv())
    return response


class Cluster(object):
    def __init__(self, seq_ids, sim_scores, taxa_separator="-", parent=None, clique=False, collapse=True, r_seed=None):
        """
        - Note that reciprocal best hits between paralogs are collapsed when instantiating group_0, so
          no problem strongly penalizing all paralogs in the scoring algorithm

        :param seq_ids: Sequence IDs
        :type seq_ids: list
        :param sim_scores: All-by-all similarity matrix for everything in the sequence_ids (and parental clusters)
        :type sim_scores: pandas.DataFrame
        :param parent: Parental sequence_ids
        :type parent: Cluster
        :param clique: Specify whether cluster is a clique or not
        :type clique: bool
        :param collapse: Specify whether reciprocal best hit paralogs should be combined
        :type collapse: bool
        :param r_seed: Set the random generator seed value
        """
        # Do an initial sanity check on the incoming graph and list of sequence ids

        expected_num_edges = int(((len(seq_ids)**2) - len(seq_ids)) / 2)
        if len(sim_scores.index) != expected_num_edges:
            raise ValueError("The number of incoming sequence ids (%s) does not match the expected graph size of %s"
                             " columns (observed %s columns)." % (len(seq_ids), expected_num_edges, len(sim_scores.index)))

        self.sim_scores = sim_scores
        self.taxa_separator = taxa_separator
        self.parent = parent

        self.subgroup_counter = 0
        # self.clique_counter = 0
        self.clique = clique
        # self.cliques = []
        self.cluster_score = None
        self.collapsed_genes = OrderedDict()  # If paralogs are reciprocal best hits, collapse them
        self.rand_gen = Random(r_seed)
        self._name = None

        self.taxa = OrderedDict()  # key = taxa id. value = list of genes coming fom that taxa.
        seq_ids = sorted(seq_ids)
        for next_seq_id in seq_ids:
            taxa = next_seq_id.split(taxa_separator)[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)

        # if clique and not parent:
        #    raise AttributeError("A clique cannot be declared without including its parental sequence_ids.")

        if parent:
            for indx, genes in parent.collapsed_genes.items():
                if indx in seq_ids:
                    self.collapsed_genes[indx] = genes
            self.seq_ids = seq_ids

        else:
            self._name = "group_0"
            # This next bit collapses all paralog reciprocal best-hit cliques so they don't gum up MCL
            # Set the full cluster score first though, in case it's needed
            self.seq_ids = seq_ids
            if collapse:
                seq_ids = self.collapse(seq_ids)

        self.seq_ids = seq_ids
        self.seq_ids_str = str(", ".join(seq_ids))
        self.seq_id_hash = helpers.md5_hash(self.seq_ids_str)

    def collapse(self, seq_ids):
        breakout = False
        while not breakout:
            breakout = True
            for seq1_id in seq_ids:
                seq1_taxa = seq1_id.split(self.taxa_separator)[0]
                paralog_best_hits = []
                for hit in self.get_best_hits(seq1_id).itertuples():
                    seq2_id = hit.seq1 if hit.seq1 != seq1_id else hit.seq2
                    if seq2_id.split(self.taxa_separator)[0] != seq1_taxa:
                        paralog_best_hits = []
                        break
                    paralog_best_hits.append(seq2_id)

                if not paralog_best_hits:
                    continue

                breakout = False
                self.collapsed_genes.setdefault(seq1_id, [])
                self.collapsed_genes[seq1_id] += paralog_best_hits
                for paralog in paralog_best_hits:
                    self.sim_scores = self.sim_scores[(self.sim_scores.seq1 != paralog) &
                                                      (self.sim_scores.seq2 != paralog)]
                    del seq_ids[seq_ids.index(paralog)]
                    del self.taxa[seq1_taxa][self.taxa[seq1_taxa].index(paralog)]
                    if paralog in self.collapsed_genes:
                        self.collapsed_genes[seq1_id] += self.collapsed_genes[paralog]
                        del self.collapsed_genes[paralog]
        return seq_ids

    def name(self):
        """
        Get cluster name, or raise error if not yet set
        :return: The cluster's name
        :rtype: str
        """
        if not self._name:
            raise AttributeError("Cluster has not been named.")
        return self._name

    def set_name(self):
        """
        Set the cluster name if not yet done, otherwise do nothing.
        The parent cluster must already have a name, and the new name here is an increment of the parents name
        :return: Nothing
        """
        if self._name:
            pass
        # elif self.clique:
        #    self._name = "%s_c%s" % (self.parent.name(), self.parent.clique_counter)
        #    self.parent.clique_counter += 1
        else:
            try:
                self.parent.name()
            except AttributeError:
                raise ValueError("Parent of current cluster has not been named.\n%s" % self)
            self._name = "%s_%s" % (self.parent.name(), self.parent.subgroup_counter)
            self.parent.subgroup_counter += 1
            # for clique in self.cliques:
            #    clique.set_name()
        return

    def compare(self, query):
        """
        Determine the percentage of sequence ids shared with a query cluster
        :param query: Another cluster object
        :return:
        """
        matches = set(self.seq_ids).intersection(query.seq_ids)
        weighted_match = (len(matches) * 2.) / (len(self) + len(query))
        print("name: %s, matches: %s, weighted_match: %s" % (self.name(), len(matches), weighted_match))
        return weighted_match

    def get_best_hits(self, gene):
        """
        Search through the sim scores data frame to find the gene(s) with the highest similarity to the query
        :param gene: Query gene name
        :type gene: str
        :return: The best hit(s). Multiple matches will all have the same similarity score.
        :rtype: pd.DataFrame
        """
        best_hits = self.sim_scores[(self.sim_scores.seq1 == gene) | (self.sim_scores.seq2 == gene)]
        if not best_hits.empty:
            best_hits = best_hits.loc[best_hits.score == max(best_hits.score)].values
            best_hits = pd.DataFrame(best_hits, columns=["seq1", "seq2", "score"])
        return best_hits

    def recursive_best_hits(self, gene, global_best_hits, tested_ids):
        """
        Compile a best-hit clique.
        The best hit for a given gene may not have the query gene as its reciprocal best hit, so find the actual best
        hit and continue until a fully contained clique is formed.
        :param gene: Query gene name
        :type gene: str
        :param global_best_hits: Growing list of best hits that is passed to each recursive call
        :type global_best_hits: pd.DataFrame
        :param tested_ids: Growing list of ids that have already been tested, so not repeating the same work
        :type tested_ids: list
        :return: The final clique
        :rtype: pd.DataFrame
        """
        best_hits = self.get_best_hits(gene)
        global_best_hits = global_best_hits.append(best_hits, ignore_index=True)
        for _edge in best_hits.itertuples():
            if _edge.seq1 not in tested_ids:
                tested_ids.append(_edge.seq1)
                global_best_hits = self.recursive_best_hits(_edge.seq1, global_best_hits, tested_ids)
            if _edge.seq2 not in tested_ids:
                tested_ids.append(_edge.seq2)
                global_best_hits = self.recursive_best_hits(_edge.seq2, global_best_hits, tested_ids)
        return global_best_hits

    def perturb(self, scores):
        """
        Add a very small bit of variance into score values when std == 0
        :param scores: Similarity scores
        :type scores: pd.DataFrame
        :return: updated data frame
        """
        valve = br.SafetyValve(global_reps=10)
        while scores.score.std() == 0:
            valve.step("Failed to perturb:\n%s" % scores)
            for indx, score in scores.score.iteritems():
                scores.set_value(indx, "score", self.rand_gen.gauss(score, (score * 0.0000001)))
        return scores

    def score(self, force=False):
        """
        :return:
        """
        if self.cluster_score and not force:
            return self.cluster_score

        unique_taxa = 0
        replicate_taxa = 0
        base_cluster = deepcopy(self.parent) if self.parent else self

        # It's possible that a cluster can have sequences not present in the parent (i.e., following orphan placement);
        # check for this, and add the sequences to base_cluster if necessary.
        items_not_in_parent = set(self.seq_ids) - set(base_cluster.seq_ids)
        for orphan in items_not_in_parent:
            base_cluster.seq_ids.append(orphan)
            orphan_taxa = orphan.split(self.taxa_separator)[0]
            base_cluster.taxa.setdefault(orphan_taxa, [])
            base_cluster.taxa[orphan_taxa].append(orphan)

        score = 1
        for taxon, genes in self.taxa.items():
            if len(genes) == 1:
                """
                When a given taxon is unique in the cluster, it gets a score relative to how many other genes are
                present in the base cluster from the same taxon. More genes == larger score.
                """
                score += len(base_cluster.taxa[taxon])
                unique_taxa += 1  # Extra improvement for larger clusters
            else:
                """
                When there is a duplicate taxon in a cluster, it gets a negative score relative to how many other genes
                are present in the base cluster form the same taxon. Fewer genes == larger negative score
                """
                replicate_taxa += len(genes) ** (2 * (len(genes) / len(base_cluster.taxa[taxon])))

        score *= (1.2 ** unique_taxa)
        if replicate_taxa:
            score /= replicate_taxa

        # Modify cluster score based on size. 'Perfect' size = half of parent
        if self.parent:  # Obviously group_0 can't be modified this way
            half_par = (len(self.parent) / 2) if len(self) <= len(self.parent) else (len(self) / 2)
            if len(self) != half_par * 2:
                score *= ((half_par - abs(len(self) - half_par)) / half_par)
            else:  # To prevent exception with log2 below
                score *= 1 / half_par
        score = log2(score)  # Convert fractions to negative scores
        score = score ** 2 if score > 0 else -1 * (score ** 2)  # Linearize the data
        self.cluster_score = score
        return self.cluster_score

    def create_rbh_cliques(self):
        results = []
        # Anything less than 4 sequences cannot be subdivided into smaller cliques
        if len(self) < 4:
            return [self]

        paralogs = OrderedDict()
        for taxa_id, group in self.taxa.items():
            if len(group) > 1:
                paralogs[taxa_id] = [self.recursive_best_hits(gene,
                                                              pd.DataFrame(columns=["seq1", "seq2", "score"]),
                                                              [gene]) for gene in group]

        # If there aren't any paralogs, we can't break up the group
        if not paralogs:
            return [self]

        # RBHCs with any overlap cannot be separated from the cluster
        seqs_to_remove = []
        for taxa_id, rbhc_dfs in paralogs.items():  # ToDo: This needs to go further, and compare any cliques between taxa_id...
            marked_for_del = []
            for i, df_i in enumerate(rbhc_dfs):
                genes_i = set(list(df_i.seq1.values) + list(df_i.seq2.values))
                if len(genes_i) < 3:  # Do not separate cliques of 2
                    marked_for_del.append(i)  # Don't need to do this for df_j, because all sets checked already
                for j, df_j in enumerate(rbhc_dfs[i+1:]):
                    genes_j = set(list(df_j.seq1.values) + list(df_j.seq2.values))
                    if list(genes_i & genes_j):
                        marked_for_del += [i, i + j]

            marked_for_del = list(set(marked_for_del))  # Remove any duplicates
            marked_for_del = sorted(marked_for_del, reverse=True)  # Delete from the largest index down
            for del_indx in marked_for_del:
                del rbhc_dfs[del_indx]
            # If all RBHCs have been disqualified, no reason to do the final test
            if not rbhc_dfs:
                continue

            # Check for overlap between similarity scores within the clique and between clique seqs and non-clique seqs
            for clique in rbhc_dfs:
                clique_ids = set(list(clique.seq1.values) + list(clique.seq2.values))
                clique_ids = list(clique_ids)
                clique_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique_ids)) &
                                                (self.sim_scores.seq2.isin(clique_ids))]
                total_scores = self.sim_scores.drop(clique_scores.index.values)

                # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
                clique_scores = self.perturb(clique_scores)
                total_scores = self.perturb(total_scores)

                total_kde = scipy.stats.gaussian_kde(total_scores.score, bw_method='silverman')
                clique_kde = scipy.stats.gaussian_kde(clique_scores.score, bw_method='silverman')
                clique_resample = clique_kde.resample(10000)
                clique95 = [np.percentile(clique_resample, 2.5), np.percentile(clique_resample, 97.5)]
                integrated = total_kde.integrate_box_1d(clique95[0], clique95[1])
                if integrated < 0.05:
                    seqs_to_remove += clique_ids
                    clique = Cluster(clique_ids, sim_scores=clique_scores, taxa_separator=self.taxa_separator,
                                     parent=self, clique=True)
                    clique.set_name()
                    results.append(clique)

        # After all that, no significant cliques to spin off as independent clusters
        if not seqs_to_remove:
            results.append(self)
            return [self]

        remaining_seqs = self.seq_ids
        for seq in set(seqs_to_remove):
            del remaining_seqs[remaining_seqs.index(seq)]

        if remaining_seqs:
            sim_scores = self.sim_scores[(self.sim_scores.seq1.isin(remaining_seqs)) &
                                         (self.sim_scores.seq2.isin(remaining_seqs))]

            remaining_cluster = Cluster(remaining_seqs, sim_scores=sim_scores, parent=self,
                                        taxa_separator=self.taxa_separator)
            remaining_cluster.set_name()
            results.append(remaining_cluster)

        return results

    """
    Much expanded scoring system that is currently broken
    def score_bak(self):
        # if self.cluster_score and not force:
        #    return self.cluster_score
        # Confirm that a score for this cluster has not been calculated before
        prev_scores = scored_clusters()
        seq_ids = helpers.md5_hash("".join(sorted(self.seq_ids)))
        # if seq_ids in prev_scores and not force:
        #    self.cluster_score = float(query(seq_ids, 'cluster_score'))
        #    return self.cluster_score

        # Don't ignore the possibility of cliques, which will alter the score.
        # Note that cliques are assumed to be the smallest unit, so not containing any sub-cliques. Valid?
        if not self.cliques:
            best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
            # pull out all genes with replicate taxa and get best hits
            for taxa, genes in self.taxa.items():
                if len(genes) > 1:
                    for gene in genes:
                        best_hits = self.recursive_best_hits(gene, best_hits, [gene])

            cliques = []
            for edge in best_hits.itertuples():
                match_indicies = []
                for indx, clique in enumerate(cliques):
                    if edge.seq1 in clique and indx not in match_indicies:
                        match_indicies.append(indx)
                    if edge.seq2 in clique and indx not in match_indicies:
                        match_indicies.append(indx)
                    if len(match_indicies) == 2:
                        break

                if not match_indicies:
                    cliques.append([edge.seq1, edge.seq2])
                elif len(match_indicies) == 1:
                    new_clique = set(cliques[match_indicies[0]] + [edge.seq1, edge.seq2])
                    cliques[match_indicies[0]] = list(new_clique)
                else:
                    match_indicies.sort()
                    new_clique = set(cliques[match_indicies[0]] + cliques[match_indicies[1]])
                    cliques[match_indicies[0]] = list(new_clique)
                    del cliques[match_indicies[1]]

            # Strip out any 'cliques' that contain less than 3 genes
            cliques = [clique for clique in cliques if len(clique) >= 3]

            # Get the similarity scores for within cliques and between cliques-remaining sequences, then generate
            # kernel-density functions for both. The overlap between the two functions is used to determine
            # whether they should be separated
            for clique in cliques:
                clique_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique)) &
                                                (self.sim_scores.seq2.isin(clique))]
                total_scores = self.sim_scores.drop(clique_scores.index.values)

                # if a clique is found that pulls in every single gene, skip
                if not len(total_scores):
                    continue

                # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
                clique_scores = self.perturb(clique_scores)
                total_scores = self.perturb(total_scores)

                total_kde = scipy.stats.gaussian_kde(total_scores.score, bw_method='silverman')
                clique_kde = scipy.stats.gaussian_kde(clique_scores.score, bw_method='silverman')
                clique_resample = clique_kde.resample(10000)
                clique95 = [np.percentile(clique_resample, 2.5), np.percentile(clique_resample, 97.5)]
                integrated = total_kde.integrate_box_1d(clique95[0], clique95[1])
                if integrated < 0.05:
                    clique = Cluster(clique, sim_scores=clique_scores, taxa_separator=self.taxa_separator,
                                     parent=self, clique=True)
                    self.cliques.append(clique)

        if self.cliques:
            clique_list = [i for j in self.cliques for i in j.seq_ids]
            decliqued_cluster = []
            for gene in self.seq_ids:
                if gene not in clique_list:
                    decliqued_cluster.append(gene)
            decliqued_cluster = decliqued_cluster
        else:
            decliqued_cluster = self.seq_ids

        self.cluster_score = self.score(decliqued_cluster)
        for clique in self.cliques:
            clique_ids = helpers.md5_hash("".join(sorted(clique.seq_ids)))
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in clique.seq_ids]))
            alb_obj = generate_msa(sb_copy)

            # if clique_ids in prev_scores:
            #    self.cluster_score += float(query(clique_ids, 'cluster_score'))
            # else:

            clique_score = self.score(clique.seq_ids)
            self.cluster_score += clique_score
            with LOCK:
                push(clique_ids, 'cluster_score', str(clique_score))
                sim_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique.seq_ids)) &
                                             (self.sim_scores.seq2.isin(clique.seq_ids))]
                push(clique_ids, 'graph', sim_scores.to_csv(index=False))
        with LOCK:
            push(seq_ids, 'cluster_score', str(self.cluster_score))
            push(seq_ids, 'graph', self.sim_scores.to_csv(index=False))
        return self.cluster_score
    """

    """The following are possible score modifier schemes to account for group size
    def gpt(self, score, taxa):  # groups per taxa
        genes_per_taxa = self.size / len(taxa.value_counts())
        num_clusters_modifier = abs(genes_per_taxa - len(self.clusters))
        score = round(score / len(self.clusters) - num_clusters_modifier ** 1.2, 2)
        return score

    def srs(self, score):  # Square root sequences
        # print(len(self.clusters))
        modifier = abs(self.size ** 0.5 - len(self.clusters)) ** 1.2
        score = round(score / len(self.clusters) - modifier, 2)
        return score
    """

    def __len__(self):
        return len(self.seq_ids)

    def __str__(self):
        return str(self.seq_ids)


def cluster2database(cluster, sql_broker, alignment):
    """
    Update the database with a cluster
    :param cluster: Cluster object
    :param sql_broker: A running helpers.SQLiteBroker object
    :param alignment: An alignment as an AlignBuddy object or string
    :return:
    """
    sql_broker.query("""INSERT OR IGNORE INTO data_table (hash, seq_ids, alignment, graph, cluster_score)
                        VALUES ('{0}', '{1}', '{2}', '{3}', '{4}')
                        """.format(cluster.seq_id_hash, cluster.seq_ids_str, str(alignment),
                                   cluster.sim_scores.to_csv(header=None, index=False), cluster.score()))
    return


# ################ PSI-PRED FUNCTIONS ################ #
def mc_psi_pred(seq_obj, args):
    outdir = args[0]
    if os.path.isfile("{0}{1}psi_pred{1}{2}.ss2".format(outdir, os.sep, seq_obj.id)):
        return
    result = run_psi_pred(seq_obj)
    with open("{0}{1}psi_pred{1}{2}.ss2".format(outdir, os.sep, seq_obj.id), "w") as ofile:
        ofile.write(result)
    return


def run_psi_pred(seq_rec):
    temp_dir = br.TempDir()
    pwd = os.getcwd()
    psipred_dir = os.path.abspath("%s%spsipred" % (os.path.dirname(__file__), os.sep))
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as ofile:
        ofile.write(seq_rec.format("fasta"))

    command = '''\
{0}{3}bin{3}seq2mtx sequence.fa > {1}{3}{2}.mtx;
{0}{3}bin{3}psipred {1}{3}{2}.mtx {0}{3}data{3}weights.dat {0}{3}data{3}weights.dat2 {0}{3}data{3}weights.dat3 > {1}{3}{2}.ss;
{0}{3}bin{3}psipass2 {0}{3}data{3}weights_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
'''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep)

    Popen(command, shell=True).wait()
    os.chdir(pwd)
    with open("%s%s%s.ss2" % (temp_dir.path, os.sep, seq_rec.id), "r") as ifile:
        result = ifile.read()
    return result


def read_ss2_file(path):
    ss_file = pd.read_csv(path, comment="#", header=None, delim_whitespace=True)
    ss_file.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
    return ss_file


def compare_psi_pred(psi1_df, psi2_df):
    num_gaps = 0
    ss_score = 0
    for row1 in psi1_df.itertuples():
        row2 = psi2_df.loc[psi2_df["indx"] == row1.indx]
        if not row2.empty:
            row_score = 0
            row_score += 1 - abs(float(row1.coil_prob) - float(row2.coil_prob))
            row_score += 1 - abs(float(row1.helix_prob) - float(row2.helix_prob))
            row_score += 1 - abs(float(row1.sheet_prob) - float(row2.sheet_prob))
            ss_score += row_score / 3
        else:
            num_gaps += 1
    align_len = len(psi2_df) + num_gaps
    ss_score /= align_len
    return ss_score
# ################ END PSI-PRED FUNCTIONS ################ #


def orthogroup_caller(master_cluster, cluster_list, seqbuddy, sql_broker, progress, outdir, psi_pred_ss2_dfs,
                      steps=1000, quiet=True, taxa_separator="-", r_seed=None):
    """
    Run MCMCMC on MCL to find the best orthogroups
    :param master_cluster: The group to be subdivided
    :type master_cluster: Cluster
    :param cluster_list: When a sequence_ids is finalized after recursion, it is appended to this list
    :param seqbuddy: The sequences that are included in the master sequence_ids
    :param sql_broker: Multithread SQL broker that can be queried
    :param progress: Progress class
    :param outdir: where are files being written to?
    :param psi_pred_ss2_dfs: OrdredDict of all ss2 dataframes with record IDs as key
    :param steps: How many MCMCMC iterations to run TODO: calculate this on the fly
    :param quiet: Suppress StdErr
    :param taxa_separator: The string that separates taxon names from gene names
    :param r_seed: Set the random generator seed value
    :return: list of sequence_ids objects
    """
    def save_cluster():
        cluster_list.append(master_cluster)
        if not os.path.isdir(os.path.join(outdir, "mcmcmc", master_cluster.name())):
            temp_dir.save(os.path.join(outdir, "mcmcmc", master_cluster.name()))
        alignment = generate_msa(seqbuddy, sql_broker)
        alignment.write(os.path.join(outdir, "alignments", master_cluster.name()))
        master_cluster.sim_scores.to_csv(os.path.join(outdir, "sim_scores", "%s.scores" % master_cluster.name()),
                                         header=None, index=False, sep="\t")
        update = len(master_cluster.seq_ids) if not master_cluster.subgroup_counter else 0
        progress.update("placed", update)
        return

    rand_gen = Random(r_seed)
    master_cluster.set_name()
    temp_dir = br.TempDir()

    # If there are no paralogs in the cluster, then it is already at its highest score and MCL is unnecessary
    keep_going = False
    for taxon, genes in master_cluster.taxa.items():
        if len(genes) > 1:
            keep_going = True
            break
    if not keep_going:
        save_cluster()
        return cluster_list
    inflation_var = mcmcmc.Variable("I", 1.1, 20, r_seed=rand_gen.randint(1, 999999999999999))
    gq_var = mcmcmc.Variable("gq", min(master_cluster.sim_scores.score), max(master_cluster.sim_scores.score),
                             r_seed=rand_gen.randint(1, 999999999999999))

    try:
        open(os.path.join(temp_dir.path, "max.txt"), "w").close()
        mcmcmc_params = ["%s" % temp_dir.path, False, seqbuddy, master_cluster,
                         taxa_separator, sql_broker, psi_pred_ss2_dfs, progress]
        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=mcmcmc_params, outfile=os.path.join(temp_dir.path, "mcmcmc_out.csv"),
                                       quiet=quiet, r_seed=rand_gen.randint(1, 999999999999999))

    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        save_cluster()
        return cluster_list

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number to start
    for chain in mcmcmc_factory.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc_factory.reset_params(["%s" % temp_dir.path, worst_score, seqbuddy, master_cluster,
                                 taxa_separator, sql_broker, psi_pred_ss2_dfs, progress])
    mcmcmc_factory.run()
    mcmcmc_output = pd.read_csv(os.path.join(temp_dir.path, "mcmcmc_out.csv"), "\t", index_col=False)
    best_score = max(mcmcmc_output["result"])
    if best_score <= master_cluster.score():
        save_cluster()
        return cluster_list

    mcl_clusters = parse_mcl_clusters(os.path.join(temp_dir.path, "best_group"))
    recursion_clusters = []
    for sub_cluster in mcl_clusters:
        cluster_ids_hash = helpers.md5_hash(", ".join(sorted(sub_cluster)))
        if len(sub_cluster) == 1:
            sim_scores = pd.DataFrame(columns=["seq1", "seq2", "score"])
        else:
            # All mcl sub clusters are written to database in mcmcmc_mcl(), so no need to check if exists
            graph = sql_broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(cluster_ids_hash))[0][0]
            sim_scores = pd.read_csv(StringIO(graph), index_col=False, header=None)
            sim_scores.columns = ["seq1", "seq2", "score"]

        sub_cluster = Cluster(sub_cluster, sim_scores=sim_scores, parent=master_cluster,
                              taxa_separator=taxa_separator, r_seed=rand_gen.randint(1, 999999999999999))
        if sub_cluster.seq_id_hash == master_cluster.seq_id_hash:
            raise ArithmeticError("The sub_cluster and master_cluster are the same, but are returning different "
                                  "scores\nsub-cluster score: %s, master score: %s\n%s"
                                  % (sub_cluster.score(), master_cluster.score(),
                                     sub_cluster.seq_id_hash))
        sub_cluster.set_name()
        if len(sub_cluster) in [1, 2]:
            cluster_list.append(sub_cluster)
            continue
        # if sub_cluster.cliques:
        #    for clique in sub_cluster.cliques:
        #        clique.clique = True
        #        clique.set_name()
        #        for _seq_id in clique.seq_ids:
        #            sub_cluster.seq_ids.remove(_seq_id)
        #    recursion_clusters += sub_cluster.cliques
        #    sub_cluster.cliques = []
        recursion_clusters.append(sub_cluster)

    for sub_cluster in recursion_clusters:
        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in sub_cluster.seq_ids])

        # Recursion... Reassign cluster_list, as all clusters are returned at the end of a call to orthogroup_caller
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, sql_broker=sql_broker,
                                         progress=progress, outdir=outdir, steps=steps, quiet=quiet,
                                         taxa_separator=taxa_separator, r_seed=rand_gen.randint(1, 999999999999999),
                                         psi_pred_ss2_dfs=psi_pred_ss2_dfs)

    save_cluster()
    return cluster_list


class Progress(object):
    def __init__(self, outdir, base_cluster):
        self.outdir = outdir
        with open(os.path.join(self.outdir, ".progress"), "w") as progress_file:
            _progress = {"mcl_runs": 0, "placed": 0, "total": len(base_cluster)}
            json.dump(_progress, progress_file)

    def update(self, key, value):
        with LOCK:
            with open(os.path.join(self.outdir, ".progress"), "r") as ifile:
                _progress = json.load(ifile)
                _progress[key] += value
            with open(os.path.join(self.outdir, ".progress"), "w") as _ofile:
                json.dump(_progress, _ofile)
        return

    def read(self):
        with LOCK:
            with open(os.path.join(self.outdir, ".progress"), "r") as ifile:
                return json.load(ifile)

    def __str__(self):
        _progress = self.read()
        return "MCL runs processed: %s. Sequences placed: %s/%s. Run time: " \
               % (_progress['mcl_runs'], _progress['placed'], _progress['total'])


def mcmcmc_mcl(args, params):
    inflation, gq, r_seed = args
    exter_tmp_dir, min_score, seqbuddy, parent_cluster, taxa_separator, sql_broker, psi_pred_ss2_dfs, progress = params
    rand_gen = Random(r_seed)
    mcl_obj = helpers.MarkovClustering(parent_cluster.sim_scores, inflation=inflation, edge_sim_threshold=gq)
    mcl_obj.run()
    progress.update('mcl_runs', 1)
    clusters = mcl_obj.clusters
    score = 0

    for indx, cluster_ids in enumerate(clusters):
        cluster_hash = helpers.md5_hash(", ".join(sorted(cluster_ids)))
        graph = sql_broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(cluster_hash))
        if graph and graph[0][0]:
            sim_scores = pd.read_csv(StringIO(graph[0][0]), index_col=False, header=None)
            score += float()
            cluster = Cluster(cluster_ids, sim_scores, parent=parent_cluster, taxa_separator=taxa_separator,
                              r_seed=rand_gen.randint(1, 999999999999999))
        else:
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in cluster_ids]))
            alb_obj = generate_msa(sb_copy, sql_broker)
            sim_scores = create_all_by_all_scores(alb_obj, psi_pred_ss2_dfs, quiet=True)
            cluster = Cluster(cluster_ids, sim_scores, parent=parent_cluster, taxa_separator=taxa_separator,
                              r_seed=rand_gen.randint(1, 999999999999999))
            cluster2database(cluster, sql_broker, alb_obj)

        clusters[indx] = cluster
        score += cluster.score()

    with LOCK:
        with open(os.path.join(exter_tmp_dir, "max.txt"), "r") as ifile:
            results = ifile.readlines()
            results = [result.strip() for result in results]
            results.append(",".join([cluster.seq_id_hash for cluster in clusters]))

        if len(results) != MCMCMC_CHAINS:
            with open(os.path.join(exter_tmp_dir, "max.txt"), "w") as ofile:
                ofile.write("\n".join(results))
        else:
            results = sorted(results)
            best_score = None
            best_clusters = []  # Hopefully just find a single best cluster, but could be more
            for clusters in results:
                score_sum = 0
                cluster_ids = []
                for cluster in clusters.split(","):
                    sql_query = sql_broker.query("SELECT seq_ids, graph FROM data_table WHERE hash='{0}'".format(cluster))
                    seq_ids = sql_query[0][0].split(", ")
                    cluster_ids.append(sql_query[0][0])
                    if len(seq_ids) > 1:  # Prevent crash if pulling a cluster with a single sequence
                        sim_scores = pd.read_csv(StringIO(sql_query[0][1]), index_col=False, header=None)
                        cluster = Cluster(seq_ids, sim_scores, parent=parent_cluster, taxa_separator=taxa_separator,
                                          r_seed=rand_gen.randint(1, 999999999999))
                        score_sum += cluster.score()

                if score_sum == best_score:
                    best_clusters.append(cluster_ids)
                elif best_score is None or score_sum > best_score:
                    best_clusters = [cluster_ids]
                    best_score = score_sum

            best_clusters = [cluster.replace(', ', '\t') for cluster in best_clusters[0]]
            with open(os.path.join(exter_tmp_dir, "best_group"), "w") as ofile:
                ofile.write('\n'.join(best_clusters))
    return score


def parse_mcl_clusters(path):
    with open(path, "r") as ifile:
        clusters = ifile.read()
    clusters = clusters.strip().split("\n")
    clusters = [cluster.strip().split("\t") for cluster in clusters]
    return clusters


def write_mcl_clusters(clusters, path):
    clusters_strings = ["\t".join(cluster.seq_ids) for cluster in clusters]
    with open(path, "w") as ofile:
        ofile.write("\n".join(clusters_strings))
    return


def check_sequences(seqbuddy, taxa_separator):
    logging.warning("Checking that the format of all sequence ids matches 'taxa%sgene'" % taxa_separator)
    failures = []
    for rec in seqbuddy.records:
        rec_id = rec.id.split(taxa_separator)
        if len(rec_id) != 2:
            failures.append(rec.id)
    if failures:
        logging.error("Malformed sequence id(s): '%s'\nThe taxa separator character is currently set to '%s',\n"
                      " which can be changed with the '-ts' flag" % (", ".join(failures), taxa_separator))
        sys.exit()
    else:
        logging.warning("    %s sequences PASSED" % len(seqbuddy))
    return


def generate_msa(seqbuddy, sql_broker):
    seq_ids = sorted([rec.id for rec in seqbuddy.records])
    seq_id_hash = helpers.md5_hash(", ".join(seq_ids))
    alignment = sql_broker.query("SELECT (alignment) FROM data_table WHERE hash='{0}'".format(seq_id_hash))
    if alignment and alignment[0][0]:
        alignment = Alb.AlignBuddy(alignment[0][0])
    else:
        if len(seqbuddy) == 1:
            alignment = Alb.AlignBuddy(str(seqbuddy))
        else:
            alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft", params="--globalpair --thread -2", quiet=True)

        sql_broker.query("""UPDATE data_table SET alignment='{0}'
                            WHERE hash='{1}'""".format(str(alignment), seq_id_hash))
    return alignment


# ################ SCORING FUNCTIONS ################ #
def create_all_by_all_scores(alignment, psi_pred_ss2_dfs, gap_open=GAP_OPEN, gap_extend=GAP_EXTEND, quiet=False):
    """
    Generate a multiple sequence alignment and pull out all-by-all similarity graph
    :param alignment: AlignBuddy object
    :param psi_pred_ss2_dfs: OrderedDict of {seqID: ss2 dataframes}
    :param gap_open: Gap initiation penalty
    :param gap_extend: Gap extension penalty
    :param quiet: Supress multicore output
    :return:
    """
    if len(alignment.records()) == 1:
        sim_scores = pd.DataFrame(data=None, columns=["seq1", "seq2", "score"])
        return sim_scores

    # Don't want to modify the alignbuddy object in place
    alignment = Alb.make_copy(alignment)

    # Need to specify what columns the PsiPred files map to now that there are gaps.
    psi_pred_ss2_dfs = deepcopy(psi_pred_ss2_dfs)  # Don't modify in place...

    for rec in alignment.records_iter():
        ss_file = psi_pred_ss2_dfs[rec.id]
        ss_counter = 0
        for indx, residue in enumerate(rec.seq):
            if residue != "-":
                psi_pred_ss2_dfs[rec.id].set_value(ss_counter, "indx", indx)
                ss_counter += 1
        psi_pred_ss2_dfs[rec.id] = ss_file

    # Scores seem to be improved by removing gaps. Need to test this explicitly for the paper though
    alignment = Alb.trimal(alignment, "gappyout")

    # Re-update PsiPred files, now that some columns are removed
    for rec in alignment.records_iter():
        new_psi_pred = [0 for _ in range(len(psi_pred_ss2_dfs[rec.id].index))]  # Instantiate list of max possible size
        indx = 0
        for row in psi_pred_ss2_dfs[rec.id].itertuples():
            if alignment.alignments[0].position_map[int(row[1])][1]:
                new_psi_pred[indx] = list(row)[1:]
                indx += 1
        new_psi_pred = new_psi_pred[:indx]
        psi_pred_ss2_dfs[rec.id] = pd.DataFrame(new_psi_pred, columns=["indx", "aa", "ss", "coil_prob",
                                                                       "helix_prob", "sheet_prob"])
    ids1 = [rec.id for rec in alignment.records_iter()]
    ids2 = [rec.id for rec in alignment.records_iter()]
    all_by_all = [0 for _ in range(int((len(ids1)**2 - len(ids1)) / 2))]
    indx = 0
    for rec1 in ids1:
        del ids2[ids2.index(rec1)]
        for rec2 in ids2:
            all_by_all[indx] = (rec1, rec2)
            indx += 1

    all_by_all_outdir = br.TempDir()
    if all_by_all:
        n = ceil(len(all_by_all) / CPUS)
        all_by_all = [all_by_all[i:i + n] for i in range(0, len(all_by_all), n)] if all_by_all else []
        score_sequences_params = [alignment, psi_pred_ss2_dfs, all_by_all_outdir.path, gap_open, gap_extend]
        br.run_multicore_function(all_by_all, mc_score_sequences, score_sequences_params,
                                  quiet=quiet)
    sim_scores_file = br.TempFile()
    sim_scores_file.write("seq1,seq2,score")
    aba_root, aba_dirs, aba_files = next(os.walk(all_by_all_outdir.path))
    for aba_file in aba_files:
        with open(os.path.join(aba_root, aba_file), "r") as ifile:
            sim_scores_file.write(ifile.read())
    sim_scores = pd.read_csv(sim_scores_file.get_handle("r"), index_col=False)
    return sim_scores


def mc_score_sequences(seq_pairs, args):
    # Calculate the best possible scores, and divide by the observed scores
    alb_obj, psi_pred_ss2_dfs, output_dir, gap_open, gap_extend = args
    file_name = helpers.md5_hash(str(seq_pairs))
    ofile = open(os.path.join(output_dir, file_name), "w")
    for seq_pair in seq_pairs:
        id1, id2 = seq_pair
        id_regex = "^%s$|^%s$" % (id1, id2)
        # Alignment comparison
        alb_copy = Alb.make_copy(alb_obj)
        Alb.pull_records(alb_copy, id_regex)
        subs_mat_score = compare_pairwise_alignment(alb_copy, gap_open, gap_extend)

        # PSI PRED comparison
        ss_score = compare_psi_pred(psi_pred_ss2_dfs[id1], psi_pred_ss2_dfs[id2])
        pair_score = (ss_score * 0.3) + (subs_mat_score * 0.7)  # Magic number weights...
        ofile.write("\n%s,%s,%s" % (id1, id2, pair_score))
    ofile.close()
    return


def compare_pairwise_alignment(alb_obj, gap_open, gap_extend):
    observed_score = 0
    seq1_best = 0
    seq2_best = 0
    seq1, seq2 = alb_obj.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(alb_obj.lengths()[0]):
        aa1 = seq1.seq[aa_pos]
        aa2 = seq2.seq[aa_pos]

        if aa1 != "-":
            seq1_best += BLOSUM62[aa1, aa1]
        if aa2 != "-":
            seq2_best += BLOSUM62[aa2, aa2]

        if aa1 == "-" or aa2 == "-":
            if prev_aa1 == "-" or prev_aa2 == "-":
                observed_score += gap_extend
            else:
                observed_score += gap_open
        else:
            observed_score += BLOSUM62[aa1, aa2]
        prev_aa1 = str(aa1)
        prev_aa2 = str(aa2)

    subs_mat_score = ((observed_score / seq1_best) + (observed_score / seq2_best)) / 2
    return subs_mat_score
# ################ END SCORING FUNCTIONS ################ #


class Orphans(object):
    def __init__(self, seqbuddy, clusters, sql_broker, psi_pred_ss2_dfs, quiet=False):
        self.seqbuddy = seqbuddy
        self.clusters = clusters
        self.sql_broker = sql_broker
        self.psi_pred_ss2_dfs = psi_pred_ss2_dfs
        self.num_orphans = 0
        self.small_clusters = OrderedDict()
        self.large_clusters = OrderedDict()
        self._separate_large_small()
        self.printer = br.DynamicPrint(quiet=quiet)
        self.tmp_file = br.TempFile()
        self.all_sim_scores = pd.Series()  # This is used for a t-test later
        parent = self.clusters[0] if self.clusters[0].parent is None else self.clusters[0].parent
        graph = self.sql_broker.query("SELECT (graph) FROM data_table "
                                      "WHERE hash='{0}'".format(parent.seq_id_hash))
        graph = pd.read_csv(StringIO(graph[0][0]), index_col=False, header=None)
        graph.columns = ["seq1", "seq2", "score"]
        self.all_sim_scores = self.all_sim_scores.append(graph.score)

    def _separate_large_small(self):
        for cluster in self.clusters:
            if len(cluster.seq_ids) > 2:
                self.large_clusters[cluster.name()] = cluster
            else:
                self.small_clusters[cluster.name()] = cluster

        if not self.small_clusters:
            logging.warning("No orphaned sequences present")
        elif not self.large_clusters:
            logging.warning("All clusters have only 1 or 2 sequences present, suggesting there are issues with your data")
            logging.info(" Are there a large number of paralogs and only a small number of taxa present?")
        elif len(self.large_clusters) == 1:
            logging.warning("Only one orthogroup with >2 sequences present, suggesting there are issues with your data")
            logging.info(" Are there a large number of taxa and only a small number of orthologs present?")
        else:
            self.num_orphans = sum([len(cluster) for group_name, cluster in self.small_clusters.items()])
            logging.warning("%s orphaned sequences present" % self.num_orphans)
        return

    def _check_orphan(self, small_cluster):
        # First prepare the data for each large cluster so they can be compared
        self.tmp_file.write("%s\n" % small_cluster.seq_ids)
        data_dict = OrderedDict()
        self.tmp_file.write("\tTotal large cluster scores\n%s\n" % self.all_sim_scores)
        for group_name, large_cluster in self.large_clusters.items():
            self.tmp_file.write("\t%s\n" % large_cluster.seq_ids)
            # Read or create an alignment containing the sequences in the small and large clusters being checked
            seq_ids = sorted(large_cluster.seq_ids + small_cluster.seq_ids)
            seqbuddy = Sb.make_copy(self.seqbuddy)
            regex = "^%s$" % "$|^".join(seq_ids)
            Sb.pull_recs(seqbuddy, regex)
            alb_obj = generate_msa(seqbuddy, self.sql_broker)

            # If a graph already exists, use it. Otherwise, make a new one for this subset of sequences.
            graph = self.sql_broker.query("SELECT (graph) FROM data_table "
                                          "WHERE hash='{0}'".format(helpers.md5_hash(", ".join(seq_ids))))
            if graph and graph[0][0]:
                sim_scores = pd.read_csv(StringIO(graph[0][0]), index_col=False, header=None)
                sim_scores.columns = ["seq1", "seq2", "score"]
            else:
                sim_scores = create_all_by_all_scores(alb_obj, self.psi_pred_ss2_dfs, quiet=True)
                cluster2database(Cluster(seq_ids, sim_scores, collapse=False), self.sql_broker, alb_obj)

            # Pull out the similarity scores between just the large and small cluster sequences
            lrg2sml_group_data = sim_scores.loc[(sim_scores['seq1'].isin(small_cluster.seq_ids)) | (sim_scores['seq2'].isin(small_cluster.seq_ids))]
            lrg2sml_group_data = lrg2sml_group_data.loc[(lrg2sml_group_data['seq1'].isin(large_cluster.seq_ids)) | (lrg2sml_group_data['seq2'].isin(large_cluster.seq_ids))]
            self.tmp_file.write("\tlrg2sml_group_data:\n%s\n" % lrg2sml_group_data)

            # We will confirm that the orphans are sufficiently similar to large group to warrant consideration using t-test
            t_test = scipy.stats.ttest_ind(lrg2sml_group_data.score, self.all_sim_scores)
            self.tmp_file.write("\t%s\n\n" % str(t_test))
            # if t_test.pvalue >= 0.05:  # Therefore, when we fail to reject the null, we consider the group.
            # Convert group_data to a numpy array so sm.stats can read it
            data_dict[group_name] = (t_test.pvalue, np.array(lrg2sml_group_data.score))

        # We only need to test the large cluster with the highest average similarity score, so find that cluster.
        averages = pd.Series()
        df = pd.DataFrame(columns=['observations', 'grouplabel'])
        for group_name, group in data_dict.items():
            averages = averages.append(pd.Series(np.mean(group[1]), index=[group_name]))
            tmp_df = pd.DataFrame(group[1], columns=['observations'])
            tmp_df['grouplabel'] = group_name
            df = df.append(tmp_df)
        max_ave = averages.argmax()

        # Confirm that the largest cluster has sufficient support (t-test value)
        if data_dict[max_ave][0] <= 0.05:
            self.tmp_file.write("No Matches\n###########################\n\n")
            return False

        # Run pairwise Tukey HSD and parse the results
        # The max_ave cluster must be significantly different from all other clusters
        result = sm.stats.multicomp.pairwise_tukeyhsd(df.observations, df.grouplabel)

        test_group_pvalues = []  # To keep track of which groups have the highest support downstream
        for line in str(result).split("\n")[4:-1]:
            line = re.sub("^ *", "", line.strip())
            line = re.sub(" +", "\t", line)
            self.tmp_file.write("%s\n" % line)
            line = line.split("\t")  # Each line --> ['group1', 'group2', 'meandiff', 'lower', 'upper', 'reject]
            if max_ave in line and 'False' in line:
                break  # Insufficient support to group the gene with max_ave group
            elif max_ave in line:
                test_group_pvalues.append(abs(float(line[2])))  # Take the abs because sign isn't meaningful
        else:
            # If the 'break' command is not encountered, the gene can be grouped with the max_ave cluster
            # Return the group name and average meandiff (allow calling code to actually do the grouping)
            self.tmp_file.write("%s\n###########################\n\n" % max_ave)
            return max_ave, data_dict[max_ave][0]
        self.tmp_file.write("No Matches\n###########################\n\n")
        return False

    def place_orphans(self):
        if not self.small_clusters:
            return
        best_cluster = OrderedDict([("small_name", None), ("large_name", None), ("meandiff", 0)])
        fostered_orphans = 0
        starting_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])
        while True:
            remaining_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])
            indx = 1
            for small_name, next_small_cluster in self.small_clusters.items():  # ToDo: Can this be multicored?
                self.printer.write("%s of %s orphans remain. Checking %s/%s orphan groups" %
                                   (remaining_orphans, starting_orphans, indx, len(self.small_clusters)))
                indx += 1
                foster_score = self._check_orphan(next_small_cluster)
                if foster_score and foster_score[1] > best_cluster["meandiff"]:
                    best_cluster = OrderedDict([("small_name", small_name), ("large_name", foster_score[0]),
                                                ("meandiff", foster_score[1])])

            if best_cluster["small_name"]:
                small_name = best_cluster["small_name"]
                large_name = best_cluster["large_name"]
                self.large_clusters[large_name].seq_ids += self.small_clusters[small_name].seq_ids

                for taxon, seq_ids in self.small_clusters[small_name].taxa.items():
                    self.large_clusters[large_name].taxa.setdefault(taxon, [])
                    self.large_clusters[large_name].taxa[taxon] += seq_ids
                logging.info("\t%s added to %s" % (" and ".join(self.small_clusters[small_name].seq_ids), large_name))
                fostered_orphans += len(self.small_clusters[small_name].seq_ids)
                del self.small_clusters[small_name]
                best_cluster = OrderedDict([("small_name", None), ("large_name", None), ("meandiff", 0)])
            else:
                break
        self.printer.clear()
        self.clusters = [cluster for group, cluster in self.small_clusters.items()] + \
                        [cluster for group, cluster in self.large_clusters.items()]

        if not fostered_orphans:
            logging.warning("No orphaned sequences were placed in orthogroups.")
        elif fostered_orphans == 1:
            logging.warning("1 orphaned sequence was placed in an orthogroups.")
        else:
            logging.warning("%s orphaned sequences have found new homes in orthogroups!" % fostered_orphans)
            logging.warning("\tNote that all of the saved alignment, mcmcmc, and sim_score\n"
                            "\tfiles will only include the original sequences.")
        return

# NOTE: There used to be a support function. Check the workscripts GitHub history


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(prog="orthogroup_caller", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("sequences", help="Location of sequence file", action="store")
    parser.add_argument("outdir", action="store", default=os.path.join(os.getcwd(), "rd-mcd"),
                        help="Where should results be written?")
    parser.add_argument("-mcs", "--mcmcmc_steps", default=1000, type=int,
                        help="Specify how deeply to sample MCL parameters")
    parser.add_argument("-sr", "--supress_recursion", action="store_true",
                        help="Stop after a single round of MCL. For testing.")
    parser.add_argument("-scc", "--supress_clique_check", action="store_true",
                        help="Do not check for or break up cliques. For testing.")
    parser.add_argument("-ssf", "--supress_singlet_folding", action="store_true",
                        help="Do not check for or merge singlets. For testing.")
    parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                        help="Don't apply the gappyout algorithm to MSAs before scoring")
    parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                        type=float, default=GAP_OPEN)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=GAP_EXTEND)
    parser.add_argument("-ts", "--taxa_separator", action="store", default="-",
                        help="Specify a string that separates taxa ids from gene names")
    parser.add_argument("-rs", "--r_seed", help="Specify a random seed for repeating a specific run", type=int)
    parser.add_argument("-f", "--force", action="store_true",
                        help="Overwrite previous run")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress all output during run (only final output is returned)")

    in_args = parser.parse_args()

    logger_obj = helpers.Logger("rdmcl.log")
    logging.info("*************************** Recursive Dynamic Markov Clustering ****************************")
    logging.warning("RD-MCL version %s\n\n%s" % (VERSION, NOTICE))
    logging.info("********************************************************************************************\n")
    logging.info("Function call: %s" % " ".join(sys.argv))
    logging.warning("\n** Environment setup **")
    if not in_args.r_seed:
        in_args.r_seed = randint(1, 999999999999)
    logging.warning("Random seed used for this run: %s" % in_args.r_seed)

    logging.info("\nWorking directory: %s" % os.getcwd())
    if not os.path.isdir(in_args.outdir):
        logging.info("Output directory: %s" % in_args.outdir)
        os.makedirs(in_args.outdir)

    if not shutil.which("mafft"):
        logging.error("The 'MAFFT' program is not detected "
                      "on your system (see http://mafft.cbrc.jp/alignment/software/).")
        sys.exit()
    mafft = Popen("mafft --version", stderr=PIPE, shell=True).communicate()[1].decode()
    logging.info("\nMAFFT version: %s" % mafft.strip())

    logging.warning("\nLaunching SQLite Daemon")

    broker = helpers.SQLiteBroker(db_file=os.path.join(in_args.outdir, "sqlite_db.sqlite"))
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()
    sequences = Sb.SeqBuddy(in_args.sequences)
    check_sequences(sequences, in_args.taxa_separator)
    seq_ids_str = ", ".join(sorted([rec.id for rec in sequences.records]))
    seq_ids_hash = helpers.md5_hash(seq_ids_str)

    # Check if job has been run already
    if broker.query("SELECT (hash) FROM data_table WHERE hash='{0}'".format(seq_ids_hash)) == seq_ids_hash:
        logging.warning("RESUME: This output directory was previous used for an identical RD-MCL run.\n"
                        "        All cached resources will be reused.")

    # Make sure all the necessary directories are present and emptied of old run files
    for _path in [os.path.join(in_args.outdir, x) for x in ["", "alignments", "mcmcmc", "sim_scores", "psi_pred"]]:
        if not os.path.isdir(_path):
            logging.info("mkdir %s" % _path)
            os.makedirs(_path)
        # Delete old 'group' files/directories
        root, dirs, files = next(os.walk(_path))
        for _file in files:
            if "group" in _file:
                os.remove(os.path.join(root, _file))
        for _dir in dirs:
            if "group" in _dir:
                shutil.rmtree(os.path.join(root, _dir))

    # Move log file into output directory
    logger_obj.move_log(os.path.join(in_args.outdir, "rdmcl.log"))

    # clique_check = True if not in_args.supress_clique_check else False
    # recursion_check = True if not in_args.supress_recursion else False

    # PSIPRED
    logging.warning("\n** PSI-Pred **")
    records_missing_ss_files = []
    records_with_ss_files = []
    for record in sequences.records:
        if os.path.isfile(os.path.join(in_args.outdir, "psi_pred", "%s.ss2" % record.id)):
            records_with_ss_files.append(record.id)
        else:
            records_missing_ss_files.append(record)
    if records_missing_ss_files and len(records_missing_ss_files) != len(sequences):
        logging.info("RESUME: PSI-Pred .ss2 files found for %s sequences:" % len(records_with_ss_files))

    if records_missing_ss_files:
        logging.warning("Executing PSI-Pred on %s sequences" % len(records_missing_ss_files))
        br.run_multicore_function(records_missing_ss_files, mc_psi_pred, [in_args.outdir])
        logging.info("\t-- finished in %s --" % TIMER.split())
        logging.info("\tfiles saved to {0}{1}psi_pred{1}".format(in_args.outdir, os.sep))
    else:
        logging.warning("RESUME: All PSI-Pred .ss2 files found in {0}{1}psi_pred{1}".format(in_args.outdir, os.sep))
    psi_pred_files = [(record.id, read_ss2_file(os.path.join(in_args.outdir, "psi_pred", "%s.ss2" % record.id)))
                      for record in sequences.records]
    psi_pred_files = OrderedDict(psi_pred_files)

    # Initial alignment
    logging.warning("\n** All-by-all graph **")
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (in_args.open_penalty, in_args.extend_penalty))

    align_data = broker.query("SELECT (alignment) FROM data_table WHERE hash='{0}'".format(seq_ids_hash))
    if align_data and align_data[0][0]:
        logging.warning("RESUME: Initial multiple sequence alignment found")
        alignbuddy = Alb.AlignBuddy(align_data[0][0])
    else:
        logging.warning("Generating initial multiple sequence alignment with MAFFT")
        alignbuddy = generate_msa(sequences, broker)
        alignbuddy.write(os.path.join(in_args.outdir, "alignments", "group_0.aln"))
        logging.info("\t-- finished in %s --\n" % TIMER.split())

    graph_data = broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(seq_ids_hash))
    if graph_data and graph_data[0][0]:
        logging.warning("RESUME: Initial all-by-all similarity graph found")
        scores_data = pd.read_csv(StringIO(graph_data[0][0]), index_col=False, header=None)
        scores_data.columns = ["seq1", "seq2", "score"]

    else:
        num_comparisons = ((len(alignbuddy.alignments[0]) ** 2) - len(alignbuddy.alignments[0])) / 2
        logging.warning("Generating initial all-by-all similarity graph (%s comparisons)" % int(num_comparisons))
        logging.info(" written to: {0}{1}sim_scores{1}complete_all_by_all.scores".format(in_args.outdir, os.sep))
        scores_data = create_all_by_all_scores(alignbuddy, psi_pred_files)
        scores_data.to_csv(os.path.join(in_args.outdir, "sim_scores", "complete_all_by_all.scores"),
                           header=None, index=False, sep="\t")
        logging.info("\t-- finished in %s --\n" % TIMER.split())

    # First push the really raw first alignment in the database, without any collapsing.
    uncollapsed_group_0 = Cluster([rec.id for rec in sequences.records], scores_data,
                                  taxa_separator=in_args.taxa_separator, collapse=False, r_seed=in_args.r_seed)
    cluster2database(uncollapsed_group_0, broker, alignbuddy)

    # Then prepare the 'real' group_0 cluster
    group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data,
                              taxa_separator=in_args.taxa_separator, r_seed=in_args.r_seed)

    cluster2database(group_0_cluster, broker, alignbuddy)

    # Base cluster score
    base_score = group_0_cluster.score()
    logging.warning("Base cluster score: %s" % round(base_score, 4))
    # Reset the score of group_0 to account for possible collapse of paralogs
    if group_0_cluster.collapsed_genes:
        logging.warning("\nReciprocal best-hit cliques of paralogs have been identified in the input sequences.")
        logging.info(" A representative sequence has been selected from each clique, and the remaining")
        logging.info(" sequences will be placed in the final clusters at the end of the run.")
        with open(os.path.join(in_args.outdir, "paralog_cliques.json"), "w") as outfile:
            json.dump(group_0_cluster.collapsed_genes, outfile)
            logging.warning(" Cliques written to: %s%sparalog_cliques.json" % (in_args.outdir, os.sep))

    # taxa_count = [x.split(in_args.taxa_separator)[0] for x in group_0_cluster.seq_ids]
    # taxa_count = pd.Series(taxa_count)
    # taxa_count = taxa_count.value_counts()

    # Ortholog caller
    logging.warning("\n** Recursive MCL **")
    final_clusters = []
    progress_tracker = Progress(in_args.outdir, group_0_cluster)

    run_time = br.RunTime(prefix=progress_tracker.__str__, _sleep=0.3, final_clear=True)
    if not in_args.quiet:
        run_time.start()
    final_clusters = orthogroup_caller(group_0_cluster, final_clusters, seqbuddy=sequences, sql_broker=broker,
                                       progress=progress_tracker, outdir=in_args.outdir, steps=in_args.mcmcmc_steps,
                                       quiet=True, taxa_separator=in_args.taxa_separator, r_seed=in_args.r_seed,
                                       psi_pred_ss2_dfs=psi_pred_files,)
    final_clusters = [cluster for cluster in final_clusters if cluster.subgroup_counter == 0]
    run_time.end()

    progress_dict = progress_tracker.read()
    logging.warning("Total MCL runs: %s" % progress_dict["mcl_runs"])
    logging.warning("\t-- finished in %s --" % TIMER.split())

    '''
    # Sort out reciprocal best hit cliques
    logging.warning("\n** Processing best-hit-cliques among clustered paralogs **")
    final_cliques = []
    for clstr in final_clusters:
        if clstr.cliques and clstr.subgroup_counter == 0:
            for next_clique in clstr.cliques:
                next_clique.set_name()
                final_cliques.append(next_clique)
                for seq_id in next_clique.seq_ids:
                    clstr.seq_ids.remove(seq_id)
            clstr.cliques = None
    final_clusters += final_cliques
    logging.warning("\t%s cliques extracted" % len(final_cliques))
    '''
    # Sort out reciprocal best hit cliques
    logging.warning("\n** Processing best-hit-cliques among clustered paralogs **")
    final_cliques = []
    for clstr in final_clusters:
        final_cliques += clstr.create_rbh_cliques()

    final_clusters = final_cliques

    # logging.warning("\t%s cliques extracted" % len(final_cliques))
    # sys.exit()
    # #################### Format the clusters and output to file #################### #
    # Fold singletons and doublets back into groups. This can't be 'resumed', because it changes the clusters
    if not in_args.supress_singlet_folding:
        logging.warning("\n** Folding orphan sequences into clusters **")
        orphans = Orphans(seqbuddy=sequences, clusters=final_clusters,
                          sql_broker=broker, psi_pred_ss2_dfs=psi_pred_files)
        orphans.place_orphans()
        orphans.tmp_file.save("temp.log")
        final_clusters = orphans.clusters
        logging.warning("\t-- finished in %s --" % TIMER.split())

    if group_0_cluster.collapsed_genes:
        logging.warning("\nPlacing collapsed paralogs into their respective clusters")
        for clust in final_clusters:
            for gene_id, paralog_list in group_0_cluster.collapsed_genes.items():
                if gene_id in clust.seq_ids:
                    clust.seq_ids += paralog_list
                    clust.taxa[gene_id.split(in_args.taxa_separator)[0]].append(gene_id)
            clust.score(force=True)

    logging.warning("Preparing final_clusters.txt")
    final_score = 0
    output = ""
    while len(final_clusters) > 0:
        _max = (0, 0)
        for ind, clust in enumerate(final_clusters):
            if len(clust.seq_ids) > _max[1]:
                _max = (ind, len(clust.seq_ids))

        ind, max_clust = _max[0], final_clusters[_max[0]]
        if max_clust.subgroup_counter == 0:  # Don't want to include parents of subgroups
            output += "%s\t%s\t" % (max_clust.name(), round(max_clust.score(), 4))
            final_score += max_clust.score()
            for seq_id in sorted(max_clust.seq_ids):
                output += "%s\t" % seq_id
            output = "%s\n" % output.strip()
        del final_clusters[ind]

    logging.warning("\t-- finished in %s --" % TIMER.split())

    logging.warning("\nTotal execution time: %s" % TIMER.total_elapsed())
    with open(os.path.join(in_args.outdir, "final_clusters.txt"), "w") as outfile:
        outfile.write(output)
        logging.warning("Final score: %s" % round(final_score, 4))
        if final_score < base_score:
            logging.warning("    This is lower than the initial base score after"
                            " the inclusion of collapsed paralogs")
        logging.warning("Clusters written to: %s%sfinal_clusters.txt" % (in_args.outdir, os.sep))

    broker.close()
