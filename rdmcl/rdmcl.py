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
import time
import argparse
import sqlite3
from io import StringIO
from subprocess import Popen, PIPE, check_output, CalledProcessError
from multiprocessing import Lock, Pipe, Process
from random import gauss, Random, randint
from math import ceil, log2
from collections import OrderedDict
from copy import deepcopy
from hashlib import md5


# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo

# My packages
try:
    import mcmcmc
    import helpers
except ImportError:
    from . import mcmcmc
    from . import helpers

from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br

# Globals
SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))
try:
    git_commit = check_output(['git', '--git-dir={0}{1}..{1}.git'.format(SCRIPT_PATH, os.sep), 'rev-parse',
                               '--short', 'HEAD'], stderr=PIPE).decode().strip()
    git_commit = " (git %s)" % git_commit if git_commit else ""
    VERSION = "1.0.4%s" % git_commit
except CalledProcessError:
    VERSION = "1.0.4"

MAFFT = shutil.which("mafft") if shutil.which("mafft") \
        else os.path.join(SCRIPT_PATH, "mafft", "bin", "mafft")
NOTICE = '''\
Public Domain Notice
--------------------
This is free software; see the LICENSE for further details.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov
--------------------
'''
LOCK = Lock()
MULTICORE_LOCK = Lock()
CPUS = br.usable_cpu_count()
TIMER = helpers.Timer()
GAP_OPEN = -5
GAP_EXTEND = 0
BLOSUM62 = helpers.make_full_mat(SeqMat(MatrixInfo.blosum62))
MCMC_CHAINS = 3
DR_BASE = 0.75
GELMAN_RUBIN = 1.1
WORKER_DB = ""
HEARTBEAT_DB = ""
WORKER_PULSE = 60
MASTER_ID = None
MASTER_PULSE = 60
PSIPREDDIR = ""


ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
               "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
for aa in ambiguous_X:
    pair = sorted((aa, "X"))
    pair = tuple(pair)
    BLOSUM62[pair] = ambiguous_X[aa]


# ToDo: Maybe remove support for taxa_sep. It's complicating my life, so just impose the '-' on users?
class Cluster(object):
    def __init__(self, seq_ids, sim_scores, taxa_sep="-", parent=None, collapse=True, r_seed=None):
        """
        - Note that reciprocal best hits between paralogs are collapsed when instantiating group_0, so
          no problem strongly penalizing all paralogs in the scoring algorithm

        :param seq_ids: Sequence IDs
        :type seq_ids: list
        :param sim_scores: All-by-all similarity matrix for everything in the sequence_ids (and parental clusters)
        :type sim_scores: pandas.DataFrame
        :param parent: Parental sequence_ids
        :type parent: Cluster
        :param collapse: Specify whether reciprocal best hit paralogs should be combined
        :type collapse: bool
        :param r_seed: Set the random generator seed value
        """
        # Do an initial sanity check on the incoming graph and list of sequence ids
        expected_num_edges = int(((len(seq_ids)**2) - len(seq_ids)) / 2)
        if len(sim_scores.index) != expected_num_edges:
            raise ValueError("The number of incoming sequence ids (%s) does not match the expected graph size of %s"
                             " rows (observed %s rows)." % (len(seq_ids), expected_num_edges, len(sim_scores.index)))

        self.sim_scores = sim_scores
        self.taxa_sep = taxa_sep
        self.parent = parent

        self.subgroup_counter = 0
        self.cluster_score = None
        self.collapsed_genes = OrderedDict()  # If paralogs are reciprocal best hits, collapse them
        self.rand_gen = Random(r_seed)
        self._name = None

        self.taxa = OrderedDict()  # key = taxa id. value = list of genes coming fom that taxa.
        seq_ids = sorted(seq_ids)
        for next_seq_id in seq_ids:
            taxa = next_seq_id.split(taxa_sep)[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)

        if parent:
            self.max_genes_in_a_taxa = parent.max_genes_in_a_taxa
            for indx, genes in parent.collapsed_genes.items():
                if indx in seq_ids:
                    self.collapsed_genes[indx] = genes
            self.seq_ids = seq_ids

        else:
            self.max_genes_in_a_taxa = max([len(self.taxa[taxa]) for taxa in self.taxa])
            self._name = "group_0"
            # This next bit collapses all paralog reciprocal best-hit cliques so they don't gum up MCL
            # Set the full cluster score first though, in case it's needed
            self.seq_ids = seq_ids
            if collapse:
                self.collapse()

        self.seq_ids_str = str(", ".join(seq_ids))
        self.seq_id_hash = helpers.md5_hash(self.seq_ids_str)

    def reset_seq_ids(self, seq_ids):
        # Note that this does NOT reset sim_scores. This needs to be updated manually
        self.seq_ids = sorted(seq_ids)
        self.seq_ids_str = str(", ".join(seq_ids))
        self.seq_id_hash = helpers.md5_hash(self.seq_ids_str)

    def collapse(self):
        breakout = False
        while not breakout:
            breakout = True
            for seq1_id in self.seq_ids:
                seq1_taxa = seq1_id.split(self.taxa_sep)[0]
                paralog_best_hits = []
                for hit in self.get_best_hits(seq1_id).itertuples():
                    seq2_id = hit.seq1 if hit.seq1 != seq1_id else hit.seq2
                    if seq2_id.split(self.taxa_sep)[0] != seq1_taxa:
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
                    del self.seq_ids[self.seq_ids.index(paralog)]
                    del self.taxa[seq1_taxa][self.taxa[seq1_taxa].index(paralog)]
                    if paralog in self.collapsed_genes:
                        self.collapsed_genes[seq1_id] += self.collapsed_genes[paralog]
                        del self.collapsed_genes[paralog]
        return

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
        else:
            try:
                self.parent.name()
            except AttributeError:
                raise ValueError("Parent of current cluster has not been named.\n%s" % self)
            self._name = "%s_%s" % (self.parent.name(), self.parent.subgroup_counter)
            self.parent.subgroup_counter += 1
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
            best_hits = best_hits.loc[best_hits.raw_score == max(best_hits.raw_score)].values
            best_hits = pd.DataFrame(best_hits, columns=list(self.sim_scores.columns.values))
        return best_hits

    def recursive_best_hits(self, gene, global_best_hits, tested_ids):
        """
        Compile a best-hit clique.
        The best hit for a given gene may not have the query gene as its reciprocal best hit, so find the actual best
        hit and continue until a fully contained clique is formed.
        Note that this does not build a COG-like clique, which should grab triangles from all included taxa (I would
        prefer to fix this)
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

    def get_base_cluster(self):
        """
        Iteratively step backwards through parents until the origin is found
        :return:
        """
        cluster = self
        while True:
            if not cluster.parent:
                return cluster
            cluster = cluster.parent

    def score(self, algorithm="dem_ret", force=False):
        if self.cluster_score and not force:
            return self.cluster_score

        assert algorithm in ["dem_ret", "drp"]

        if algorithm == "dem_ret":
            return self._score_diminishing_returns()
        elif algorithm == "drp":
            return self._score_direct_replicate_penalty()

    def _score_diminishing_returns(self):
        """
        Arrange all genes into a table with species as column headings, filling in rows as possible.
        The top row will therefore have the most sequences in it, with equal or fewer in each subsequent row.
        The first row gets a full score. The second gets 1/2 score, third 1/4, nth 1/2^(n-1).
        :return:
        """
        base_cluster = deepcopy(self.get_base_cluster())

        # It's possible that a cluster can have sequences not present in the parent (i.e., following orphan placement);
        # check for this, and add the sequences to base_cluster if necessary.
        items_not_in_parent = set(self.seq_ids) - set(base_cluster.seq_ids)
        for orphan in items_not_in_parent:
            base_cluster.seq_ids.append(orphan)
            orphan_taxa = orphan.split(self.taxa_sep)[0]
            base_cluster.taxa.setdefault(orphan_taxa, [])
            base_cluster.taxa[orphan_taxa].append(orphan)

        # Split cluster up into subclusters containing no more than one gene from each taxa
        # eg. ["Mle1", "Baf1", "Lla1", "Mle2", "Pdb1", "Lla2", "Mle3"] would become:
        #     [["Mle1", "Baf1", "Lla1", "Pdb1"], ["Mle2", "Lla2"], ["Mle3"]]
        subclusters = [[]]
        for taxon, genes in self.taxa.items():
            for indx, gene in enumerate(genes):
                if len(subclusters) > indx:
                    subclusters[indx].append(taxon)
                else:
                    subclusters.append([taxon])

        score = 0
        for indx, subcluster in enumerate(subclusters):
            subscore = 0
            # Minimum score for a single gene is 1, which is given to the genes contained in the the taxa with the
            # largest number of genes. All other genes are pegged to this value, and increase in value if included in
            # a cluster proportionally to how many genes are in the taxa.
            for taxon in subcluster:
                subscore += self.max_genes_in_a_taxa / len(base_cluster.taxa[taxon])

            # Multiply subscore by 1-2, based on how many taxa contained (perfect score if all possible taxa present)
            subscore *= len(subcluster) / len(base_cluster.taxa) + 1

            # Reduce subscores as we encounter replicate taxa
            subscore *= DR_BASE ** indx

            score += subscore

        self.cluster_score = score
        return self.cluster_score

    def _score_direct_replicate_penalty(self):
        # This is currently only accessible by modifying the default in score()
        unique_taxa = 0
        replicate_taxa = 0
        base_cluster = deepcopy(self.parent) if self.parent else self

        # It's possible that a cluster can have sequences not present in the parent (i.e., following orphan placement);
        # check for this, and add the sequences to base_cluster if necessary.
        items_not_in_parent = set(self.seq_ids) - set(base_cluster.seq_ids)
        for orphan in items_not_in_parent:
            base_cluster.seq_ids.append(orphan)
            orphan_taxa = orphan.split(self.taxa_sep)[0]
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
            else:  # To prevent exception with log2 below when parent and child are same length
                score *= 1 / half_par
        score = log2(score)  # Convert fractions to negative scores
        score = score ** 2 if score > 0 else -1 * (score ** 2)  # Linearize the data
        self.cluster_score = score
        return self.cluster_score

    def pull_scores_subgraph(self, seq_ids):
        assert not set(seq_ids) - set(self.seq_ids)  # Confirm that all requested seq_ids are present
        subgraph = self.sim_scores[(self.sim_scores.seq1.isin(seq_ids)) & (self.sim_scores.seq2.isin(seq_ids))]
        return subgraph

    def create_rbh_cliques(self, log_file=None):
        """
        Check for reciprocal best hit cliques
        :param log_file: Handle to a writable file
        :return: list of Cluster objects or False
        """
        log_file = log_file if log_file else br.TempFile()
        log_file.write("# ####### Testing %s ####### #\n%s\n" % (self.name(), self.seq_ids))
        results = []
        # Anything less than 6 sequences cannot be subdivided into smaller cliques, because it would create orphans
        if len(self) < 6:
            log_file.write("\tTERMINATED: Group too small.\n\n")
            return [self]

        paralogs = OrderedDict()
        for taxa_id, group in self.taxa.items():
            if len(group) > 1:
                for gene in group:
                    rbhc = self.recursive_best_hits(gene,
                                                    pd.DataFrame(columns=list(self.sim_scores.columns.values)),
                                                    [gene])
                    # If ANY clique encompasses the entire cluster, we're done
                    if len(set(list(rbhc.seq1.values) + list(rbhc.seq2.values))) == len(self):
                        log_file.write("\tTERMINATED: Entire cluster pulled into clique on %s." % taxa_id)
                        return [self]

                    paralogs.setdefault(taxa_id, [])
                    paralogs[taxa_id].append(rbhc)

        # If there aren't any paralogs, we can't break up the group
        if not paralogs:
            log_file.write("\tTERMINATED: No paralogs present.\n\n")
            return [self]

        # RBHCs with any overlap within a taxon cannot be separated from the cluster
        seqs_to_remove = []
        log_file.write("\tChecking taxa for overlapping cliques:\n")
        for taxa_id, rbhc_dfs in paralogs.items():
            log_file.write("\t\t# #### %s #### #\n" % taxa_id)
            marked_for_del = OrderedDict()
            for i, df_i in enumerate(rbhc_dfs):
                genes_i = set(list(df_i.seq1.values) + list(df_i.seq2.values))
                # Do not separate cliques of 2. df_j sizes are covered in df_i loop.
                if len(genes_i) < 3:
                    marked_for_del[i] = "is too small"
                for j, df_j in enumerate(rbhc_dfs[i+1:]):
                    genes_j = set(list(df_j.seq1.values) + list(df_j.seq2.values))
                    if list(genes_i & genes_j):
                        if i not in marked_for_del:
                            marked_for_del[i] = "overlaps with %s" % sorted(list(genes_j))
                        if i + j not in marked_for_del:
                            marked_for_del[i + j] = "overlaps with %s" % sorted(list(genes_i))

            if marked_for_del:
                log_file.write("\t\tDisqualified cliques:\n")
                # Re-order del indices to delete from the largest index down
                marked_for_del = sorted(marked_for_del.items(), key=lambda x: x[0], reverse=True)
                marked_for_del = OrderedDict(marked_for_del)
                for del_indx, reason in marked_for_del.items():
                    clique_ids = set(list(rbhc_dfs[del_indx].seq1.values) + list(rbhc_dfs[del_indx].seq2.values))
                    log_file.write("\t\t\t%s %s\n" % (sorted(list(clique_ids)), reason))
                    del rbhc_dfs[del_indx]

            # If all RBHCs have been disqualified, no reason to do the final test
            if not rbhc_dfs:
                log_file.write("\n\t\t!! ALL CLIQUES DISQUALIFIED !!\n\n")
                continue

            # Check for overlap between similarity scores within the clique and between clique seqs and non-clique seqs
            log_file.write("\t\tTest for KDE separation:\n")
            for clique in rbhc_dfs:
                clique_ids = set(list(clique.seq1.values) + list(clique.seq2.values))
                clique_ids = list(clique_ids)
                log_file.write("\t\t\t%s\n" % sorted(clique_ids))
                clique_scores = self.pull_scores_subgraph(clique_ids)
                outer_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique_ids)) |
                                               (self.sim_scores.seq2.isin(clique_ids))]

                outer_scores = outer_scores[((outer_scores.seq1.isin(clique_ids))
                                             & (outer_scores.seq2.isin(clique_ids) == False)) |  # Note the == must stay
                                            ((outer_scores.seq2.isin(clique_ids))
                                             & (outer_scores.seq1.isin(clique_ids) == False))]  # Note the == must stay

                # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
                clique_scores = self.perturb(clique_scores)
                outer_scores = self.perturb(outer_scores)

                total_kde = scipy.stats.gaussian_kde(outer_scores.raw_score, bw_method='silverman')
                log_file.write("""\
\t\t\tOuter KDE: {'shape': %s, 'covariance': %s, 'inv_cov': %s, '_norm_factor': %s}
""" % (total_kde.dataset.shape, total_kde.covariance[0][0], total_kde.inv_cov[0][0], total_kde._norm_factor))

                clique_kde = scipy.stats.gaussian_kde(clique_scores.raw_score, bw_method='silverman')
                log_file.write("""\
\t\t\tClique KDE: {'shape': %s, 'covariance': %s, 'inv_cov': %s,  '_norm_factor': %s}
""" % (clique_kde.dataset.shape, clique_kde.covariance[0][0], clique_kde.inv_cov[0][0], clique_kde._norm_factor))

                clique_resample = clique_kde.resample(10000)  # ToDo: figure out how to control this with r_seed!
                clique95 = [np.percentile(clique_resample, 2.5), np.percentile(clique_resample, 97.5)]
                log_file.write("\t\t\tclique95: %s\n" % clique95)

                integrated = total_kde.integrate_box_1d(clique95[0], clique95[1])
                log_file.write("\t\t\tintegrated: %s\n" % integrated)
                if integrated < 0.05:
                    log_file.write("\t\t\tPASS\n\n")
                    seqs_to_remove += clique_ids
                    results.append([clique_ids, clique_scores])
                else:
                    log_file.write("\t\t\tFAIL\n\n")

        # After all that, no significant cliques to spin off as independent clusters...
        if not seqs_to_remove:
            log_file.write("\tTERMINATED: No significant cliques identified.\n\n")
            return [self]

        # Look for any overlap between the cliques identified for each taxon and merge them together
        if len(results) > 1:
            log_file.write("\t\tChecking if new cliques can combine\n")
        combined = []
        while results:
            clique_i = results[0]
            drop_j_indx = []
            for indx, clique_j in enumerate(results[1:]):
                log_file.write("\t\t%s - %s  " % (sorted(clique_i[0]), sorted(clique_j[0])))
                if list(set(clique_i[0]) & set(clique_j[0])):
                    log_file.write("YES\n")
                    clique_i[0] = set(clique_i[0] + clique_j[0])
                    clique_i[0] = sorted(list(clique_i[0]))
                    clique_i[1] = self.pull_scores_subgraph(clique_i[0])
                    drop_j_indx.append(indx + 1)
                else:
                    log_file.write("NO\n")
            combined.append(clique_i)
            for indx in sorted(drop_j_indx, reverse=True):
                del results[indx]
            del results[0]
        remaining_seqs = self.seq_ids
        for seq in set(seqs_to_remove):
            del remaining_seqs[remaining_seqs.index(seq)]

        # Do not leave orphans by spinning off cliques
        if len(remaining_seqs) in [1, 2]:
            log_file.write("\tTERMINATED: Spinning off cliques would orphan %s.\n\n" % " and ".join(remaining_seqs))
            return [self]

        for indx, result in enumerate(combined):
            cluster_ids, cluster_sim_scores = result
            cluster = Cluster(cluster_ids, sim_scores=cluster_sim_scores, taxa_sep=self.taxa_sep,
                              parent=self)
            cluster.set_name()
            results.append(cluster)

        if remaining_seqs:
            sim_scores = self.pull_scores_subgraph(remaining_seqs)

            remaining_cluster = Cluster(remaining_seqs, sim_scores=sim_scores, parent=self,
                                        taxa_sep=self.taxa_sep)
            remaining_cluster.set_name()
            results.append(remaining_cluster)
        log_file.write("\tCliques identified and spun off:\n\t\t%s\n\n" %
                       "\n\t\t".join([str(res.seq_ids) for res in results]))
        return results

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
    if os.path.isfile("{0}{1}{2}.ss2".format(outdir, os.sep, seq_obj.id)):
        return
    result = run_psi_pred(seq_obj)
    with open("{0}{1}{2}.ss2".format(outdir, os.sep, seq_obj.id), "w") as ofile:
        ofile.write(result)
    return


def run_psi_pred(seq_rec):
    temp_dir = br.TempDir()
    pwd = os.getcwd()
    psipred_dir = os.path.abspath("%s%spsipred" % (SCRIPT_PATH, os.sep))
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as ofile:
        ofile.write(seq_rec.format("fasta"))

    if shutil.which("psipred"):
        command = '''\
seq2mtx sequence.fa > {1}{3}{2}.mtx;
psipred {1}{3}{2}.mtx {0}{3}data{3}weights.dat {0}{3}data{3}weights.dat2 {0}{3}data{3}weights.dat3 > {1}{3}{2}.ss;
psipass2 {0}{3}data{3}weights_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
'''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep)

    else:
        data_weights = "{0}{1}data{1}weights".format(psipred_dir, os.sep)
        command = '''\
    {0}{3}bin{3}seq2mtx sequence.fa > {1}{3}{2}.mtx;
    {0}{3}bin{3}psipred {1}{3}{2}.mtx {4}.dat {4}.dat2 {4}.dat3 > {1}{3}{2}.ss;
    {0}{3}bin{3}psipass2 {4}_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
    '''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep, data_weights)

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
                      steps=1000, chains=3, walkers=2, quiet=True, taxa_sep="-", r_seed=None, convergence=None):
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
    :param chains: Number of MCMCMC chains to spin off
    :param walkers: Number of Metropolis-Hastings walkers per chain
    :param quiet: Suppress StdErr
    :param taxa_sep: The string that separates taxon names from gene names
    :param r_seed: Set the random generator seed value
    :param convergence: Set minimum Gelman-Rubin PSRF value for convergence
    :return: list of sequence_ids objects
    """
    def save_cluster(end_message=None):
        cluster_list.append(master_cluster)
        if not os.path.isfile(os.path.join(temp_dir.path, "best_group")):
            with open(os.path.join(temp_dir.path, "best_group"), "w") as _ofile:
                _ofile.write('\t'.join(master_cluster.seq_ids))
        if end_message:
            with open(os.path.join(temp_dir.path, "end_message.log"), "w") as _ofile:
                _ofile.write(end_message + "\n")

        if not os.path.isdir(os.path.join(outdir, "mcmcmc", master_cluster.name())):
            temp_dir.save(os.path.join(outdir, "mcmcmc", master_cluster.name()))
        _, alignment = create_all_by_all_scores(seqbuddy, psi_pred_ss2_dfs, sql_broker, quiet=True)
        alignment.write(os.path.join(outdir, "alignments", master_cluster.name()))
        master_cluster.sim_scores.to_csv(os.path.join(outdir, "sim_scores", "%s.scores" % master_cluster.name()),
                                         header=None, index=False, sep="\t")
        update = len(master_cluster.seq_ids) if not master_cluster.subgroup_counter else 0
        progress.update("placed", update)
        return

    rand_gen = Random(r_seed)
    master_cluster.set_name()
    temp_dir = br.TempDir()
    convergence = GELMAN_RUBIN if convergence is None else float(convergence)

    # If there are no paralogs in the cluster, then it is already at its highest score and MCL is unnecessary
    keep_going = False
    for taxon, genes in master_cluster.taxa.items():
        if len(genes) > 1:
            keep_going = True
            break
    if not keep_going:
        save_cluster("No paralogs")
        return cluster_list
    inflation_var = mcmcmc.Variable("I", 1.1, 20, r_seed=rand_gen.randint(1, 999999999999999))
    gq_var = mcmcmc.Variable("gq", min(master_cluster.sim_scores.score), max(master_cluster.sim_scores.score),
                             r_seed=rand_gen.randint(1, 999999999999999))

    try:
        open(os.path.join(temp_dir.path, "max.txt"), "w").close()
        mcmcmc_params = ["%s" % temp_dir.path, seqbuddy, master_cluster,
                         taxa_sep, sql_broker, psi_pred_ss2_dfs, progress]
        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1, quiet=quiet,
                                       num_walkers=walkers, num_chains=chains, convergence=convergence,
                                       outfiles=os.path.join(temp_dir.path, "mcmcmc_out"), params=mcmcmc_params,
                                       include_lava=True, include_ice=True, r_seed=rand_gen.randint(1, 999999999999999))

    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        save_cluster("MCMCMC failed to find parameters")
        return cluster_list

    # I know what the best and worst possible scores are, so let MCMCMC know (better for calculating acceptance rates)
    # The worst score possible would be all genes in each taxa segregated.
    worst_score = 0
    for taxon, seq_ids in master_cluster.taxa.items():
        bad_clust = Cluster(seq_ids, master_cluster.pull_scores_subgraph(seq_ids), parent=master_cluster)
        worst_score += bad_clust.score()

    # The best score would be perfect separation of taxa
    sub_clusters = [[]]
    for taxon, genes in master_cluster.taxa.items():
        for indx, gene in enumerate(genes):
            if len(sub_clusters) > indx:
                sub_clusters[indx].append(gene)
            else:
                sub_clusters.append([gene])
    best_score = 0
    for seq_ids in sub_clusters:
        subcluster = Cluster(seq_ids, master_cluster.pull_scores_subgraph(seq_ids), parent=master_cluster)
        best_score += subcluster.score()

    for chain in mcmcmc_factory.chains:
        for walker in chain.walkers:
            walker.score_history += [worst_score, best_score]

    mcmcmc_factory.reset_params(["%s" % temp_dir.path, seqbuddy, master_cluster,
                                 taxa_sep, sql_broker, psi_pred_ss2_dfs, progress])
    mcmcmc_factory.run()
    best_score = pd.DataFrame()
    for indx in range(len(mcmcmc_factory.chains)):
        mcmcmc_output = pd.read_csv(os.path.join(temp_dir.path, "mcmcmc_out_%s.csv" % (indx+1)), "\t", index_col=False)
        result = mcmcmc_output.loc[mcmcmc_output["result"] == mcmcmc_output["result"].max()]
        best_score = result if best_score.empty or result["result"].iloc[0] > best_score["result"].iloc[0] \
            else best_score

    if round(best_score["result"].iloc[0], 8) <= round(master_cluster.score(), 8):
        save_cluster("New best score of %s is less than master cluster at %s"
                     % (best_score["result"].iloc[0], master_cluster.score()))
        return cluster_list

    mcl_obj = helpers.MarkovClustering(master_cluster.sim_scores, inflation=best_score["I"].iloc[0],
                                       edge_sim_threshold=best_score["gq"].iloc[0])
    mcl_obj.run()
    progress.update('mcl_runs', 1)
    mcl_clusters = mcl_obj.clusters

    # Write out the actual best clusters
    best_clusters = ['\t'.join(cluster) for cluster in mcl_clusters]
    with open(os.path.join(temp_dir.path, "best_group"), "w") as ofile:
        ofile.write('\n'.join(best_clusters))

    recursion_clusters = []
    for sub_cluster in mcl_clusters:
        cluster_ids_hash = helpers.md5_hash(", ".join(sorted(sub_cluster)))
        if len(sub_cluster) == 1:
            sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
        else:
            # All mcl sub clusters are written to database in mcmcmc_mcl(), so no need to check if exists
            graph = sql_broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(cluster_ids_hash))[0][0]
            sim_scores = pd.read_csv(StringIO(graph), index_col=False, header=None)
            sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]

        sub_cluster = Cluster(sub_cluster, sim_scores=sim_scores, parent=master_cluster,
                              taxa_sep=taxa_sep, r_seed=rand_gen.randint(1, 999999999999999))
        if sub_cluster.seq_id_hash == master_cluster.seq_id_hash:  # This shouldn't ever happen
            raise ArithmeticError("The sub_cluster and master_cluster are the same, but are returning different "
                                  "scores\nsub-cluster score: %s, master score: %s\n%s"
                                  % (best_score["result"].iloc[0], master_cluster.score(),
                                     sub_cluster.seq_id_hash))
        sub_cluster.set_name()
        if len(sub_cluster) in [1, 2]:
            cluster_list.append(sub_cluster)
            continue
        recursion_clusters.append(sub_cluster)

    for sub_cluster in recursion_clusters:
        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in sub_cluster.seq_ids])

        sub_cluster.collapse()

        # Recursion... Reassign cluster_list, as all clusters are returned at the end of a call to orthogroup_caller
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, sql_broker=sql_broker,
                                         progress=progress, outdir=outdir, steps=steps, quiet=quiet, chains=chains,
                                         walkers=walkers, taxa_sep=taxa_sep, convergence=convergence,
                                         r_seed=rand_gen.randint(1, 999999999999999), psi_pred_ss2_dfs=psi_pred_ss2_dfs)

    save_cluster("Sub clusters returned")
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
    """
    Function passed to mcmcmcm.MCMCMC that will execute MCL and return the final cluster scores
    :param args: Sample values to run MCL with and a random seed [inflation, gq, r_seed]
    :param params: List of parameters (see below for unpacking assignment)
    :return:
    """
    inflation, gq, r_seed = args
    exter_tmp_dir, seqbuddy, parent_cluster, taxa_sep, sql_broker, psi_pred_ss2_dfs, progress = params
    rand_gen = Random(r_seed)
    mcl_obj = helpers.MarkovClustering(parent_cluster.sim_scores, inflation=inflation, edge_sim_threshold=gq)
    mcl_obj.run()
    progress.update('mcl_runs', 1)
    clusters = mcl_obj.clusters
    score = 0

    for indx, cluster_ids in enumerate(clusters):
        sb_copy = Sb.make_copy(seqbuddy)
        sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in cluster_ids]))
        sim_scores, alb_obj = create_all_by_all_scores(sb_copy, psi_pred_ss2_dfs, sql_broker, quiet=True)
        cluster = Cluster(cluster_ids, sim_scores, parent=parent_cluster, taxa_sep=taxa_sep,
                          r_seed=rand_gen.randint(1, 999999999999999))

        clusters[indx] = cluster
        score += cluster.score()

    with LOCK:
        with open(os.path.join(exter_tmp_dir, "max.txt"), "r") as ifile:
            results = ifile.readlines()
            results = [result.strip() for result in results]
            results.append(",".join([cluster.seq_id_hash for cluster in clusters]))
            results = sorted(results)

        with open(os.path.join(exter_tmp_dir, "max.txt"), "w") as ofile:
            ofile.write("\n".join(results))

        if len(results) == MCMC_CHAINS:
            best_score = None
            best_clusters = []  # Hopefully just find a single best set of cluster, but could be more
            for clusters in results:
                score_sum = 0
                cluster_ids = []
                for cluster in clusters.split(","):
                    sql_query = sql_broker.query("SELECT seq_ids, graph FROM data_table WHERE hash='%s'" % cluster)
                    seq_ids = sql_query[0][0].split(", ")
                    cluster_ids.append(sql_query[0][0])
                    if len(seq_ids) > 1:  # Prevent crash if pulling a cluster with a single sequence
                        sim_scores = pd.read_csv(StringIO(sql_query[0][1]), index_col=False, header=None)
                        cluster = Cluster(seq_ids, sim_scores, parent=parent_cluster, taxa_sep=taxa_sep,
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


def check_sequences(seqbuddy, taxa_sep):
    logging.warning("Checking that the format of all sequence ids matches 'taxa%sgene'" % taxa_sep)
    failures = []
    for rec in seqbuddy.records:
        rec_id = rec.id.split(taxa_sep)
        if len(rec_id) != 2:
            failures.append(rec.id)
    if failures:
        logging.error("Malformed sequence id(s): '%s'\nThe taxa separator character is currently set to '%s',\n"
                      " which can be changed with the '-ts' flag" % (", ".join(failures), taxa_sep))
        sys.exit()
    else:
        logging.warning("    %s sequences PASSED" % len(seqbuddy))
    return


class HeartBeat(object):
    def __init__(self, hbdb_path, pulse_rate):
        self.hbdb_path = hbdb_path
        self.pulse_rate = pulse_rate
        self.master_id = None
        self.running_process = None
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            try:
                cursor.execute('CREATE TABLE heartbeat (thread_id INTEGER PRIMARY KEY AUTOINCREMENT, '
                               'thread_type TEXT, pulse INTEGER)')
            except sqlite3.OperationalError:
                pass

    def _run(self, check_file_path):
        split_time = time.time()
        while True:
            with open("%s" % check_file_path, "r") as ifile:
                ifile_content = ifile.read()
            if ifile_content != "Running":
                break

            if split_time < time.time() - self.pulse_rate:
                with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
                    cursor.execute("UPDATE heartbeat SET pulse=%s WHERE thread_id=%s" % (round(time.time()),
                                                                                         self.master_id))
                split_time = time.time()
            time.sleep(randint(1, 100) / 100)
        return

    def start(self):
        if self.running_process:
            self.end()
        tmp_file = br.TempFile()
        tmp_file.write("Running")
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            cursor.execute("INSERT INTO heartbeat (thread_type, pulse) "
                           "VALUES ('master', %s)" % round(time.time()))
            self.master_id = cursor.lastrowid
        p = Process(target=self._run, args=(tmp_file.path,))
        p.daemon = 1
        p.start()
        self.running_process = [tmp_file, p]
        return

    def end(self):
        if not self.running_process:
            return
        self.running_process[0].clear()
        while self.running_process[1].is_alive():
            continue
        self.running_process = None
        with helpers.ExclusiveConnect(self.hbdb_path) as cursor:
            cursor.execute("DELETE FROM heartbeat WHERE thread_id=%s" % self.master_id)
        self.master_id = None
        return


# ################ SCORING FUNCTIONS ################ #
def create_all_by_all_scores(seqbuddy, psi_pred_ss2_dfs, sql_broker, gap_open=GAP_OPEN,
                             gap_extend=GAP_EXTEND, quiet=False):
    """
    Generate a multiple sequence alignment and pull out all-by-all similarity graph
    :param seqbuddy: AlignBuddy object
    :param psi_pred_ss2_dfs: OrderedDict of {seqID: ss2 dataframes}
    :param sql_broker: Active broker object to search/update SQL database
    :param gap_open: Gap initiation penalty
    :param gap_extend: Gap extension penalty
    :param quiet: Supress multicore output
    :return:
    """
    seq_ids = sorted([rec.id for rec in seqbuddy.records])
    seq_id_hash = helpers.md5_hash(", ".join(seq_ids))

    if len(seqbuddy) == 1:
        sim_scores = pd.DataFrame(data=None, columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
        alignment = Alb.AlignBuddy(str(seqbuddy))
        cluster2database(Cluster(seq_ids, sim_scores, collapse=False), sql_broker, alignment)
        return sim_scores, alignment

    # Grab from the database first, if the data exists there already
    query = sql_broker.query("SELECT graph, alignment FROM data_table WHERE hash='%s'" % seq_id_hash)
    if query and len(query[0]) == 2:
        sim_scores, alignment = query[0]
        sim_scores = pd.read_csv(StringIO(sim_scores), index_col=False, header=None)
        sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
        return sim_scores, Alb.AlignBuddy(alignment)

    # Try to feed the job to independent workers
    if WORKER_DB and os.path.isfile(WORKER_DB):
        heartbeat = HeartBeat(HEARTBEAT_DB, MASTER_PULSE)
        heartbeat.start()

        run_worker = True
        already_queued = False
        # complete_check = []
        queue_check = []
        processing_check = []
        heartbeat_check = []

        with helpers.ExclusiveConnect(HEARTBEAT_DB) as cursor:
            workers = cursor.execute("SELECT * FROM heartbeat "
                                     "WHERE thread_type='worker' "
                                     "AND pulse>%s" % (time.time() - WORKER_PULSE)).fetchone()

        if workers:
            with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                # Make sure the job hasn't already been submitted or completed)
                complete_check = cursor.execute("SELECT hash FROM complete WHERE hash='%s'" % seq_id_hash).fetchone()
                queue_check = cursor.execute("SELECT hash FROM queue WHERE hash='%s'" % seq_id_hash).fetchone()
                processing_check = cursor.execute("SELECT hash FROM processing WHERE hash='%s'"
                                                  % seq_id_hash).fetchone()

                if not complete_check and not queue_check and not processing_check:
                    cursor.execute("INSERT INTO queue (hash, seqs, psi_pred_dir, master_id) "
                                   "VALUES ('%s', '%s', '%s', %s)" %
                                   (seq_id_hash, str(seqbuddy), PSIPREDDIR, heartbeat.master_id))
                cursor.execute("INSERT INTO waiting (hash, master_id) VALUES ('%s', %s)"
                               % (seq_id_hash, heartbeat.master_id))
        else:
            run_worker = False

        vanished = 0
        while run_worker:
            query = sql_broker.query("SELECT graph, alignment FROM data_table WHERE hash='%s'" % seq_id_hash)
            if query and len(query[0]) == 2:
                sim_scores, alignment = query[0]
                sim_scores = pd.read_csv(StringIO(sim_scores), index_col=False, header=None)
                sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
                with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                    waiting_check = cursor.execute("SELECT * FROM waiting WHERE hash='%s'" % seq_id_hash).fetchall()
                    if len(waiting_check) == 1:
                        cursor.execute("DELETE FROM complete WHERE hash='%s'" % seq_id_hash)
                    cursor.execute("DELETE FROM waiting WHERE hash='%s' AND master_id=%s"
                                   % (seq_id_hash, heartbeat.master_id))
                heartbeat.end()
                return sim_scores, Alb.AlignBuddy(alignment)

            ######################

            with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                complete_check = cursor.execute("SELECT * FROM complete WHERE hash='%s'" % seq_id_hash).fetchone()
                if not complete_check:
                    # Grab data and release database lock
                    queue_check = cursor.execute("SELECT * FROM queue "
                                                 "WHERE hash='%s'" % seq_id_hash).fetchone()
                    processing_check = cursor.execute("SELECT * FROM processing "
                                                      "WHERE hash='%s'" % seq_id_hash).fetchone()

                else:
                    waiting_check = cursor.execute("SELECT * FROM waiting WHERE hash='%s'" % seq_id_hash).fetchall()
                    if len(waiting_check) == 1:
                        cursor.execute("DELETE FROM complete WHERE hash='%s'" % seq_id_hash)
                    cursor.execute("DELETE FROM waiting WHERE hash='%s' AND master_id=%s"
                                   % (seq_id_hash, heartbeat.master_id))
                    alignment = Alb.AlignBuddy(StringIO(complete_check[1]))
                    sim_scores = pd.read_csv(StringIO(complete_check[2]), index_col=False, header=None)
                    sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
                    cluster2database(Cluster(seq_ids, sim_scores, collapse=False), sql_broker, alignment)
                    heartbeat.end()
                    return sim_scores, alignment

            with helpers.ExclusiveConnect(HEARTBEAT_DB) as cursor:
                heartbeat_check = cursor.execute("SELECT * FROM heartbeat "
                                                 "WHERE thread_type='worker' "
                                                 "AND pulse>%s" % (time.time() - (WORKER_PULSE * 3))).fetchall()

            # Ensure there are still workers kicking around
            if not heartbeat_check:
                run_worker = False
                with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                    cursor.execute("DELETE FROM queue")
                    cursor.execute("DELETE FROM processing")
                    cursor.execute("DELETE FROM waiting WHERE master_id=%s" % heartbeat.master_id)
            else:
                if processing_check:
                    # Make sure the respective worker is still alive to finish processing this job
                    workers = [(thread_id, pulse) for thread_id, thread_type, pulse in heartbeat_check]
                    workers = OrderedDict(workers)
                    worker_id = processing_check[1]
                    if worker_id not in workers:
                        # Doesn't look like it... Might as well queue it back up
                        with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                            try:
                                cursor.execute("INSERT INTO queue (hash, seqs, psi_pred_dir, master_id) "
                                               "VALUES ('%s', '%s', '%s', %s)" %
                                               (seq_id_hash, str(seqbuddy), PSIPREDDIR, heartbeat.master_id))
                                cursor.execute("DELETE FROM processing WHERE hash='%s'" % seq_id_hash)
                            except sqlite3.IntegrityError as err:
                                # Ooops, someone else must have beat us too it
                                if "UNIQUE constraint failed: queue.hash" not in str(err):
                                    raise err
                elif queue_check:
                    pass
                else:
                    # Ummmm... Where'd the job go??? Queue it back up if it really seems to have vanished
                    if vanished == 5:
                        print("We had a vanish!")
                        vanished = 0
                        with helpers.ExclusiveConnect(WORKER_DB) as cursor:
                            try:
                                cursor.execute("INSERT INTO queue (hash, seqs, psi_pred_dir, master_id) "
                                               "VALUES ('%s', '%s', '%s', %s)" %
                                               (seq_id_hash, str(seqbuddy), PSIPREDDIR, heartbeat.master_id))
                            except sqlite3.IntegrityError as err:
                                # Ooops, someone else must have beat us too it
                                if "UNIQUE constraint failed: queue.hash" not in str(err):
                                    raise err
                    else:
                        vanished += 1
                time.sleep(randint(1, 100) / 100)  # Pause for some part of one second
        heartbeat.end()

    # If the job couldn't be pushed off on a worker, do it directly
    alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), MAFFT, params="--globalpair --thread -1", quiet=True)

    # Need to specify what columns the PsiPred files map to now that there are gaps.
    psi_pred_ss2_dfs = deepcopy(psi_pred_ss2_dfs)  # Don't modify in place...

    for rec in seqbuddy.records:
        ss_file = psi_pred_ss2_dfs[rec.id]
        ss_counter = 0
        for indx, residue in enumerate(rec.seq):
            if residue != "-":
                psi_pred_ss2_dfs[rec.id].set_value(ss_counter, "indx", indx)
                ss_counter += 1
        psi_pred_ss2_dfs[rec.id] = ss_file

    # Scores seem to be improved by removing gaps. ToDo: Need to test this explicitly for the paper
    alignment = Alb.trimal(alignment, threshold=0.9)

    # Re-update PsiPred files now that some columns, possibly including non-gap characters, are removed
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
        n = int(ceil(len(all_by_all) / CPUS))
        all_by_all = [all_by_all[i:i + n] for i in range(0, len(all_by_all), n)] if all_by_all else []
        score_sequences_params = [alignment, psi_pred_ss2_dfs, all_by_all_outdir.path, gap_open, gap_extend]
        with MULTICORE_LOCK:
            br.run_multicore_function(all_by_all, mc_score_sequences, score_sequences_params,
                                      quiet=quiet, max_processes=CPUS)
    sim_scores_file = br.TempFile()
    sim_scores_file.write("seq1,seq2,subsmat,psi")
    aba_root, aba_dirs, aba_files = next(os.walk(all_by_all_outdir.path))
    for aba_file in aba_files:
        with open(os.path.join(aba_root, aba_file), "r") as ifile:
            sim_scores_file.write(ifile.read())
    sim_scores = pd.read_csv(sim_scores_file.get_handle("r"), index_col=False)

    # Set raw score, which is used by Orphan placement
    sim_scores['raw_score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
    # Distribute final scores_components between 0-1.
    sim_scores['psi'] = (sim_scores['psi'] - sim_scores['psi'].min()) / \
                        (sim_scores['psi'].max() - sim_scores['psi'].min())

    sim_scores['subsmat'] = (sim_scores['subsmat'] - sim_scores['subsmat'].min()) / \
                            (sim_scores['subsmat'].max() - sim_scores['subsmat'].min())

    # ToDo: Experiment testing these magic number weights...
    sim_scores['score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
    cluster2database(Cluster(seq_ids, sim_scores, collapse=False), sql_broker, alignment)
    return sim_scores, alignment


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
        alb_copy = Sb.SeqBuddy(alb_copy.records())
        Sb.clean_seq(alb_copy)
        # assert len(psi_pred_ss2_dfs[id1]) == len(alb_copy.records[0].seq)
        # assert len(psi_pred_ss2_dfs[id2]) == len(alb_copy.records[1].seq)
        ss_score = compare_psi_pred(psi_pred_ss2_dfs[id1], psi_pred_ss2_dfs[id2])
        ofile.write("\n%s,%s,%s,%s" % (id1, id2, subs_mat_score, ss_score))
    ofile.close()
    return


def compare_pairwise_alignment(alb_obj, gap_open, gap_extend):
    observed_score = 0
    observed_len = alb_obj.lengths()[0]
    seq1_best = 0
    seq1_len = 0
    seq2_best = 0
    seq2_len = 0
    seq1, seq2 = alb_obj.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(observed_len):
        aa1 = seq1.seq[aa_pos]
        aa2 = seq2.seq[aa_pos]

        if aa1 != "-":
            seq1_best += BLOSUM62[aa1, aa1]
            seq1_len += 1
        if aa2 != "-":
            seq2_best += BLOSUM62[aa2, aa2]
            seq2_len += 1
        if aa1 == "-" or aa2 == "-":
            if prev_aa1 == "-" or prev_aa2 == "-":
                observed_score += gap_extend
            else:
                observed_score += gap_open
        else:
            observed_score += BLOSUM62[aa1, aa2]
        prev_aa1 = str(aa1)
        prev_aa2 = str(aa2)

    # Calculate average per-residue log-odds ratios for both best possible alignments and observed
    # Note: The best score range is 4 to 11. Possible observed range is -4 to 11.
    observed_score = (observed_score / observed_len)

    seq1_best = (seq1_best / seq1_len)
    seq1_score = observed_score / seq1_best

    seq2_best = (seq2_best / seq2_len)
    seq2_score = observed_score / seq2_best

    subs_mat_score = (seq1_score + seq2_score) / 2
    return subs_mat_score
# ################ END SCORING FUNCTIONS ################ #


class Orphans(object):
    def __init__(self, seqbuddy, clusters, sql_broker, psi_pred_ss2_dfs, quiet=False):
        """
        Organizes all of the orphan prediction/folding logic into a class
        :param seqbuddy: This should include all of the sequences in the entire population being tests
        :param clusters: List of all cluster objects the sequences have currently been grouped into
        :param sql_broker: Open SQLiteBroker object
        :param psi_pred_ss2_dfs: OrderedDict of all PSI Pred graphs (seq_id: df)
        :param quiet: Suppress any terminal output
        """
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

        # Create a null model from all within cluster scores from the large clusters; used for placement tests later
        self.lrg_cluster_sim_scores = pd.Series()
        for clust_name, clust in self.large_clusters.items():
            self.lrg_cluster_sim_scores = self.lrg_cluster_sim_scores.append(clust.sim_scores.raw_score,
                                                                             ignore_index=True)

    def _separate_large_small(self):
        for cluster in self.clusters:
            if len(cluster.seq_ids) > 2:
                self.large_clusters[cluster.name()] = cluster
            else:
                self.small_clusters[cluster.name()] = cluster

    def temp_merge_clusters(self, clust1, clust2):
        # Read or create an alignment containing the sequences in the clusters being merged
        seq_ids = sorted(clust1.seq_ids + clust2.seq_ids)
        seqbuddy = Sb.make_copy(self.seqbuddy)
        regex = "^%s$" % "$|^".join(seq_ids)
        Sb.pull_recs(seqbuddy, regex)

        sim_scores, alb_obj = create_all_by_all_scores(seqbuddy, self.psi_pred_ss2_dfs, self.sql_broker, quiet=True)
        cluster2database(Cluster(seq_ids, sim_scores, collapse=False), self.sql_broker, alb_obj)
        return sim_scores

    def _check_orphan(self, small_cluster):
        # First prepare the data for each large cluster so they can be compared
        log_output = "%s\n" % small_cluster.seq_ids
        log_output += "Min sim score needed: %s\n" % self.lrg_cluster_sim_scores.min()
        data_dict = OrderedDict()
        for group_name, large_cluster in self.large_clusters.items():
            log_output += "\t%s\n" % group_name
            log_output += "\t%s\n" % large_cluster.seq_ids
            log_output += "\tAmong large:\t%s - %s\n" % (round(large_cluster.sim_scores.raw_score.min(), 4),
                                                         round(large_cluster.sim_scores.raw_score.max(), 4))

            sim_scores = self.temp_merge_clusters(large_cluster, small_cluster)

            # Pull out the similarity scores between just the large and small cluster sequences
            lrg2sml_group_data = sim_scores.loc[(sim_scores['seq1'].isin(small_cluster.seq_ids))
                                                | (sim_scores['seq2'].isin(small_cluster.seq_ids))]
            lrg2sml_group_data = lrg2sml_group_data.loc[(lrg2sml_group_data['seq1'].isin(large_cluster.seq_ids))
                                                        | (lrg2sml_group_data['seq2'].isin(large_cluster.seq_ids))]
            log_output += "\tLarge 2 small:\t%s - %s\n\n" % (round(lrg2sml_group_data.raw_score.min(), 4),
                                                             round(lrg2sml_group_data.raw_score.max(), 4))

            # Confirm that the orphans are sufficiently similar to large group to warrant consideration by ensuring
            # that at least one similarity score is greater than the lowest score with all large groups.
            # Also, convert group_data to a numpy array so sm.stats can read it
            if lrg2sml_group_data.raw_score.max() >= self.lrg_cluster_sim_scores.min():
                data_dict[group_name] = (True, np.array(lrg2sml_group_data.raw_score))
            else:
                data_dict[group_name] = (False, np.array(lrg2sml_group_data.raw_score))
        # We only need to test the large cluster with the highest average similarity score, so find that cluster.
        averages = pd.Series()
        df = pd.DataFrame(columns=['observations', 'grouplabel'])
        for group_name, group in data_dict.items():
            averages = averages.append(pd.Series(np.mean(group[1]), index=[group_name]))
            tmp_df = pd.DataFrame(group[1], columns=['observations'])
            tmp_df['grouplabel'] = group_name
            df = df.append(tmp_df)
        max_ave_name = averages.argmax()

        # Confirm that the largest cluster has sufficient support
        if not data_dict[max_ave_name][0]:
            log_output += "No Matches: Best cluster (%s) insufficient support\n###########################\n\n"\
                          % max_ave_name
            with LOCK:
                self.tmp_file.write(log_output)
            return False

        groupsunique, groupintlab = np.unique(df.grouplabel, return_inverse=True)
        if len(groupsunique) < 2:
            # The gene can be grouped with the max_ave cluster because it's the only large cluster available
            log_output += "\n%s added to %s\n###########################\n\n" % (small_cluster.seq_ids, max_ave_name)
            mean_diff = abs(np.mean(data_dict[max_ave_name][1]) - np.mean(self.lrg_cluster_sim_scores))
            return max_ave_name, mean_diff

        # Run pairwise Tukey HSD and parse the results
        # The max_ave cluster must be significantly different from all other clusters
        result = sm.stats.multicomp.pairwise_tukeyhsd(df.observations, df.grouplabel)

        success = True
        for line in str(result).split("\n")[4:-1]:
            line = re.sub("^ *", "", line.strip())
            line = re.sub(" +", "\t", line)
            line = line.split("\t")  # Each line --> ['group1', 'group2', 'meandiff', 'lower', 'upper', 'reject]
            if max_ave_name in line:
                log_output += "%s\n" % "\t".join(line)
                if 'False' in line:
                    success = False  # Insufficient support to group the gene with max_ave group
        if success:
            # The gene can be grouped with the max_ave cluster
            # Return the group name and average meandiff (allow calling code to actually do the grouping)
            log_output += "\n%s added to %s\n###########################\n\n" % (small_cluster.seq_ids, max_ave_name)
            mean_diff = abs(np.mean(data_dict[max_ave_name][1]) - np.mean(self.lrg_cluster_sim_scores))
            with LOCK:
                self.tmp_file.write(log_output)
            return max_ave_name, mean_diff

        log_output += "Best group (%s) fails Tukey HSD\n###########################\n\n" % max_ave_name
        with LOCK:
            self.tmp_file.write(log_output)
        return False

    def mc_check_orphans(self, small_cluster, args):
        tmp_file = args[0]
        foster_score = self._check_orphan(small_cluster) if self.large_clusters else False
        if foster_score:
            large_name, mean_diff = foster_score
            with LOCK:
                with open(tmp_file, "a") as ofile:
                    ofile.write("%s\t%s\t%s\n" % (small_cluster.name(), large_name, mean_diff))
        return

    def place_orphans(self, multi_core=True):
        if not self.small_clusters:
            return
        best_cluster = OrderedDict([("small_name", None), ("large_name", None), ("meandiff", 0)])
        fostered_orphans = 0
        starting_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])

        # Loop through repeatedly, adding a single small cluster at a time and updating the large clusters on each loop
        while True and self.large_clusters:
            remaining_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])
            self.printer.write("%s of %s orphans remain" % (remaining_orphans, starting_orphans))
            tmp_file = br.TempFile()
            small_clusters = [clust for clust_name, clust in self.small_clusters.items()]

            # Create each large/small combination and add to database if necessary. This has overhead here, but prevents
            # over subscribing available CPU resources in mc_check_orphans()
            for sm_clust in small_clusters:
                for _, lg_clust in self.large_clusters.items():
                    self.temp_merge_clusters(sm_clust, lg_clust)

            if multi_core:
                with MULTICORE_LOCK:
                    br.run_multicore_function(small_clusters, self.mc_check_orphans, [tmp_file.path],
                                              quiet=True, max_processes=CPUS)
            else:
                for clust in small_clusters:
                    self.mc_check_orphans(clust, [tmp_file.path])

            lines = tmp_file.read()
            if lines:
                lines = lines.strip().split("\n")
                lines = sorted(lines)
                for line in lines:
                    small_name, large_name, mean_diff = line.split("\t")
                    mean_diff = float(mean_diff)
                    if mean_diff > best_cluster["meandiff"]:
                        best_cluster = OrderedDict([("small_name", small_name), ("large_name", large_name),
                                                    ("meandiff", mean_diff)])

            if best_cluster["small_name"]:
                small_name = best_cluster["small_name"]
                large_name = best_cluster["large_name"]
                large_cluster = self.large_clusters[large_name]
                small_cluster = self.small_clusters[small_name]
                large_cluster.reset_seq_ids(small_cluster.seq_ids + large_cluster.seq_ids)

                seqs = Sb.pull_recs(Sb.make_copy(self.seqbuddy), "^%s$" % "$|^".join(large_cluster.seq_ids))
                sim_scores, alb_obj = create_all_by_all_scores(seqs, self.psi_pred_ss2_dfs,
                                                               self.sql_broker, quiet=True)

                large_cluster.sim_scores = sim_scores
                for taxon, seq_ids in small_cluster.taxa.items():
                    large_cluster.taxa.setdefault(taxon, [])
                    large_cluster.taxa[taxon] += seq_ids
                    for seq_id, paralogs in small_cluster.collapsed_genes.items():
                        large_cluster.collapsed_genes[seq_id] = paralogs
                fostered_orphans += len(small_cluster.seq_ids)
                del self.small_clusters[small_name]
                del small_cluster
                best_cluster = OrderedDict([("small_name", None), ("large_name", None), ("meandiff", 0)])
            else:
                break

        self.printer.clear()
        self.clusters = [cluster for group, cluster in self.small_clusters.items()] + \
                        [cluster for group, cluster in self.large_clusters.items()]

        return


def argparse_init():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="rdmcl", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                     description='''\
\033[1mRD-MCL\033[m
  Got orthogroups?

  Recursive Dynamic Markov Clustering identifies MCL parameters that
  maximize  one-to-one  placement  of  sequences  into  orthogroups.
  Differences in evolutionary rates between orthogroups is accounted
  for  by  recursively  decomposing  large clusters  where possible.

  \033[1m\033[91mAll  sequences  must  be  homologous! RD-MCL  is \033[4mNOT\033[24m  intended for
  identifying  orthogroups  from  whole  genome/transcriptome  data.\033[m

\033[1mUsage\033[m:
  rdmcl "/path/to/sequence_file" [outdir] [-options]
''')

    parser.register('action', 'setup', _SetupAction)

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("sequences", help="Location of sequence file (most formats are fine)", action="store")
    outdir = os.path.join(os.getcwd(), "rdmcd-%s" % time.strftime("%d-%m-%Y"))
    positional.add_argument("outdir", action="store", help="Where should results be written? (default=%s)" % outdir,
                            nargs="?", default=outdir)

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("-rs", "--r_seed", type=int, metavar="",
                              help="Specify a random seed for repeating a specific run")
    parser_flags.add_argument("-sql", "--sqlite_db", action="store", metavar="",
                              help="If RD-MCL has already created a SQLite database, reusing it can speed things up")
    parser_flags.add_argument("-psi", "--psipred_dir", action="store", metavar="",
                              help="If RD-MCL has already calculated PSI-Pred files, point to the directory")
    parser_flags.add_argument("-ts", "--taxa_sep", action="store", default="-", metavar="",
                              help="Specify the string that separates taxa ids from gene names (default='-')")
    parser_flags.add_argument("-ch", "--chains", default=MCMC_CHAINS, type=int, metavar="",
                              help="Specify how many MCMC chains to run (default=3)")
    parser_flags.add_argument("-wlk", "--walkers", default=3, type=int, metavar="",
                              help="Specify how many Metropolis-Hastings walkers are in each chain (default=3)")
    parser_flags.add_argument("-drb", "--dr_base", type=float, metavar="",
                              help="Set the base for the diminishing-returns "
                                   "orthogroup scoring function (default=%s)" % DR_BASE)
    parser_flags.add_argument("-cnv", "--converge", type=float, metavar="",
                              help="Set minimum Gelman-Rubin PSRF value for convergence (default=%s)" % GELMAN_RUBIN)
    parser_flags.add_argument("-cpu", "--max_cpus", type=int, action="store", default=CPUS, metavar="",
                              help="Specify the maximum number of cores RD-MCL can use (default=%s)" % CPUS)
    parser_flags.add_argument("-lwt", "--lock_wait_time", type=int, default=1200, metavar="",
                              help="Specify num seconds a process should wait on the SQLite database before crashing"
                                   " out (default=1200)")
    parser_flags.add_argument("-mcs", "--mcmc_steps", default=0, type=int, metavar="",
                              help="Specify a max number of MCMC steps (default=auto-detect)")
    parser_flags.add_argument("-op", "--open_penalty", type=float, default=GAP_OPEN, metavar="",
                              help="Penalty for opening a gap in pairwise alignment scoring (default=%s)" % GAP_OPEN)
    parser_flags.add_argument("-ep", "--ext_penalty", type=float, default=GAP_EXTEND, metavar="",
                              help="Penalty to extend a gap in pairwise alignment scoring (default=%s)" % GAP_EXTEND)
    parser_flags.add_argument("-wdb", "--workdb", action="store", default="",
                              help="Specify the directory that independent workers will be monitoring")
    parser_flags.add_argument("-f", "--force", action="store_true",
                              help="Try to run no matter what.")
    parser_flags.add_argument("-q", "--quiet", action="store_true",
                              help="Suppress output during run (only final output is returned)")

    # Developer testing
    dev_flags = parser.add_argument_group(title="\033[1mDeveloper commands (Caution!)\033[m")
    dev_flags.add_argument("-spc", "--suppress_paralog_collapse", action="store_true",
                           help="Do not merge best hit paralogs")
    dev_flags.add_argument("-sr", "--suppress_recursion", action="store_true",
                           help="Stop after a single round of MCL")
    dev_flags.add_argument("-scc", "--suppress_clique_check", action="store_true",
                           help="Do not check for or break up cliques")
    dev_flags.add_argument("-ssf", "--suppress_singlet_folding", action="store_true",
                           help="Do not check for or merge singlets")
    dev_flags.add_argument("-sit", "--suppress_iteration", action="store_true",
                           help="Only check for cliques and orphans once")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")
    misc.add_argument('-v', '--version', action='version', version="RD-MCL version %s\n\n%s" % (VERSION, NOTICE))
    misc.add_argument("-setup", action="setup", dest=argparse.SUPPRESS, default=argparse.SUPPRESS)

    in_args = parser.parse_args()
    return in_args


# Create a new argparse action type to allow the setup option to override the positional arguments
class _SetupAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 _help="Complete the RDMCL installation"):
        super(_SetupAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=_help)

    def __call__(self, parser, namespace, values, option_string=None):
        setup()
        parser.exit()


def full_run(in_args):
    logger_obj = helpers.Logger("rdmcl.log")
    logging.info("*************************** Recursive Dynamic Markov Clustering ****************************")
    logging.warning("RD-MCL version %s\n\n%s" % (VERSION, NOTICE))
    logging.info("********************************************************************************************\n")

    if not os.path.isfile(os.path.join(SCRIPT_PATH, "config.ini")):
        print("""Error: You have not run the rdmcl setup script.
Please do so now:

    $: %s -setup

""" % __file__)
        return

    logging.info("Function call: %s" % " ".join(sys.argv))
    logging.info("Date: %s" % time.strftime("%H:%M-%d-%m-%Y"))
    logging.warning("\n** Environment setup **")
    if not in_args.r_seed:
        in_args.r_seed = randint(1, 999999999999)
    logging.warning("Random seed used for this run: %s\n" % in_args.r_seed)

    logging.info("Working directory: %s" % os.getcwd())
    logging.warning("Output directory:\n    %s" % os.path.abspath(in_args.outdir))
    if not os.path.isdir(in_args.outdir):
        os.makedirs(in_args.outdir)

    if not os.path.isfile(MAFFT):
        logging.error("The 'MAFFT' program is not detected "
                      "on your system. Please run \n\t$: %s -setup" % __file__)
        sys.exit()

    mafft = Popen("%s --version" % MAFFT, stderr=PIPE, shell=True).communicate()[1].decode()
    logging.info("\nMAFFT version: %s" % mafft.strip())
    logging.info("SeqBuddy version: %s.%s" % (Sb.VERSION.major, Sb.VERSION.minor))
    logging.info("AlignBuddy version: %s.%s" % (Alb.VERSION.major, Sb.VERSION.minor))

    logging.warning("\nLaunching SQLite Daemons")

    sqlite_path = os.path.join(in_args.outdir, "sqlite_db.sqlite") if not in_args.sqlite_db \
        else os.path.abspath(in_args.sqlite_db)
    broker = helpers.SQLiteBroker(db_file=sqlite_path, lock_wait_time=in_args.lock_wait_time)
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()
    sequences = Sb.SeqBuddy(in_args.sequences)

    global WORKER_DB
    global HEARTBEAT_DB

    if in_args.workdb:
        in_args.workdb = os.path.abspath(in_args.workdb)
        WORKER_DB = os.path.join(in_args.workdb, "work_db.sqlite")

        HEARTBEAT_DB = os.path.join(in_args.workdb, "heartbeat_db.sqlite")
        heartbeat = HeartBeat(HEARTBEAT_DB, MASTER_PULSE)
        heartbeat.start()

    # Do a check for large data sets and confirm run
    if len(sequences) > 1000 and not in_args.force:
        msg = """\
WARNING: You are attempting to run RD-MCL on %s sequences. This will take a long time!
Be aware that RD-MCL is NOT intended to identify orthogroups from among non-homologous
sequences. If this is your goal, you should use other software like OrthoMCL or ProteinOrtho.
RD-MCL is specifically for fine grained classification of individual protein families.
Continue? y/[n] """ % len(sequences)
        logging.info("\n%s\n" % msg)

        if not br.ask(msg, default="no", timeout=120):
            logging.warning("\nRun terminated. If you really do want to run this large job, use the -f flag.")
            sys.exit()
        else:
            logging.warning("Proceeding with large run.\n")

    sequences = Sb.delete_metadata(sequences)
    check_sequences(sequences, in_args.taxa_sep)
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

    # Prepare log files into output directory
    logger_obj.move_log(os.path.join(in_args.outdir, "rdmcl.log"))
    if os.path.isfile(os.path.join(in_args.outdir, "orphans.log")):
        os.remove(os.path.join(in_args.outdir, "orphans.log"))

    if os.path.isfile(os.path.join(in_args.outdir, "cliques.log")):
        os.remove(os.path.join(in_args.outdir, "cliques.log"))

    with open(os.path.join(in_args.outdir, "paralog_cliques"), "w") as outfile:
        outfile.write("###########################################################\n"
                      "# If a named cluster contains reciprocal best hit cliques #\n"
                      "# among a group of paralogs, they are collapsed down to a #\n"
                      "# single representative. The collapses are itemized here. #\n"
                      "###########################################################\n\n")

    # Set CPU limits
    global CPUS
    CPUS = in_args.max_cpus

    # PSIPRED
    logging.warning("\n** PSI-Pred **")
    records_missing_ss_files = []
    records_with_ss_files = []
    if in_args.psipred_dir and not os.path.isdir(in_args.psipred_dir):
        logging.warning("PSI-Pred directory indicated below does not exits:\n%s \tUsing default location instead."
                        % in_args.psipred_dir)
        in_args.psipred_dir = False

    in_args.psipred_dir = "{0}{1}psi_pred".format(in_args.outdir, os.sep) if not in_args.psipred_dir \
        else in_args.psipred_dir

    global PSIPREDDIR
    PSIPREDDIR = os.path.abspath(in_args.psipred_dir)

    for record in sequences.records:
        if os.path.isfile(os.path.join(in_args.psipred_dir, "%s.ss2" % record.id)):
            records_with_ss_files.append(record.id)
        else:
            records_missing_ss_files.append(record)
    if records_missing_ss_files and len(records_missing_ss_files) != len(sequences):
        logging.info("RESUME: PSI-Pred .ss2 files found for %s sequences:" % len(records_with_ss_files))

    if records_missing_ss_files:
        logging.warning("Executing PSI-Pred on %s sequences" % len(records_missing_ss_files))
        br.run_multicore_function(records_missing_ss_files, mc_psi_pred, [in_args.psipred_dir], max_processes=CPUS)
        logging.info("\t-- finished in %s --" % TIMER.split())
        logging.info("\tfiles saved to {0}{1}".format(in_args.psipred_dir, os.sep))
    else:
        logging.warning("RESUME: All PSI-Pred .ss2 files found")

    psi_pred_files = []
    for record in sequences.records:
        psi_pred_files.append((record.id, read_ss2_file(os.path.join(in_args.psipred_dir, "%s.ss2" % record.id))))

    psi_pred_files = OrderedDict(psi_pred_files)

    # Initial alignment
    logging.warning("\n** All-by-all graph **")
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (in_args.open_penalty, in_args.ext_penalty))

    """
    align_data = broker.query("SELECT (alignment) FROM data_table WHERE hash='{0}'".format(seq_ids_hash))
    if align_data and align_data[0][0]:
        logging.warning("RESUME: Initial multiple sequence alignment found")
        alignbuddy = Alb.AlignBuddy(align_data[0][0])
    else:
        logging.warning("Generating initial multiple sequence alignment with MAFFT")
        alignbuddy = generate_msa(sequences, broker)
        alignbuddy.write(os.path.join(in_args.outdir, "alignments", "group_0.aln"))
        logging.info("\t-- finished in %s --\n" % TIMER.split())
    """

    num_comparisons = ((len(sequences) ** 2) - len(sequences)) / 2
    logging.warning("Generating initial all-by-all similarity graph (%s comparisons)" % int(num_comparisons))
    logging.info(" written to: {0}{1}sim_scores{1}complete_all_by_all.scores".format(in_args.outdir, os.sep))
    scores_data, alignbuddy = create_all_by_all_scores(sequences, psi_pred_files, broker)
    scores_data.to_csv(os.path.join(in_args.outdir, "sim_scores", "complete_all_by_all.scores"),
                       header=None, index=False, sep="\t")
    logging.info("\t-- finished in %s --\n" % TIMER.split())

    if in_args.dr_base:
        global DR_BASE
        DR_BASE = in_args.dr_base
    logging.warning("Diminishing returns scoring base: %s" % DR_BASE)

    # First prepare the really raw first alignment in the database, without any collapsing.
    uncollapsed_group_0 = Cluster([rec.id for rec in sequences.records], scores_data,
                                  taxa_sep=in_args.taxa_sep, collapse=False, r_seed=in_args.r_seed)

    # Then prepare the 'real' group_0 cluster
    if not in_args.suppress_paralog_collapse:
        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data,
                                  taxa_sep=in_args.taxa_sep, r_seed=in_args.r_seed)
    else:
        group_0_cluster = uncollapsed_group_0

    cluster2database(group_0_cluster, broker, alignbuddy)

    # Base cluster score
    base_score = group_0_cluster.score()
    logging.warning("Base cluster score: %s" % round(base_score, 4))

    # Ortholog caller
    logging.warning("\n** Recursive MCL **")
    global GELMAN_RUBIN
    if in_args.mcmc_steps >= 100:
        logging.warning("User specified maximum MCMCMC steps: %s" % in_args.mcmc_steps)
    elif in_args.mcmc_steps == 0:
        logging.warning("Auto detect MCMCMC convergence")
    else:
        logging.warning("User specified maximum MCMCMC steps of %s is too low (min=100). "
                        "Switching to auto-detect" % in_args.mcmc_steps)
        in_args.mcmc_steps = 0

    if in_args.mcmc_steps < 100:
        GELMAN_RUBIN = in_args.converge if in_args.converge else GELMAN_RUBIN
        logging.warning("Gelman-Rubin convergence breakpoint: %s" % GELMAN_RUBIN)

    if in_args.chains >= 2:
        logging.warning("Number of MCMC chains: %s" % in_args.chains)
    else:
        logging.warning("User specified value of %s MCMC chains is too low. "
                        "Switching to 3" % in_args.chains)
        in_args.chains = 3

    global MCMC_CHAINS
    MCMC_CHAINS = in_args.chains

    if in_args.walkers >= 2:
        logging.warning("Number of Metropolis-Hastings walkers per chain: %s" % in_args.walkers)
    else:
        logging.warning("User specified value of %s Metropolis-Hastings walkers per chain is too low. "
                        "Switching to 2" % in_args.walkers)
        in_args.walkers = 2

    final_clusters = []
    progress_tracker = Progress(in_args.outdir, group_0_cluster)

    run_time = br.RunTime(prefix=progress_tracker.__str__, _sleep=0.3, final_clear=True)
    if not in_args.quiet:
        run_time.start()
    final_clusters = orthogroup_caller(group_0_cluster, final_clusters, seqbuddy=sequences, sql_broker=broker,
                                       progress=progress_tracker, outdir=in_args.outdir, steps=in_args.mcmc_steps,
                                       quiet=True, taxa_sep=in_args.taxa_sep, r_seed=in_args.r_seed,
                                       psi_pred_ss2_dfs=psi_pred_files, chains=MCMC_CHAINS, walkers=in_args.walkers,
                                       convergence=GELMAN_RUBIN)
    final_clusters = [cluster for cluster in final_clusters if cluster.subgroup_counter == 0]
    run_time.end()

    progress_dict = progress_tracker.read()
    logging.warning("Total MCL runs: %s" % progress_dict["mcl_runs"])
    logging.warning("\t-- finished in %s --" % TIMER.split())

    if not in_args.suppress_clique_check or not in_args.suppress_singlet_folding:
        logging.warning("\n** Iterative placement of orphans and paralog RBHC removal **")
    else:
        logging.warning("\n** Skipping placement of orphans and paralog RBHC removal **")
    check_md5 = None
    while check_md5 != md5(str(sorted([clust.seq_ids for clust in final_clusters])).encode("utf-8")).hexdigest():
        check_md5 = md5(str(sorted([clust.seq_ids for clust in final_clusters])).encode("utf-8")).hexdigest()

        # Sort out reciprocal best hit cliques
        if not in_args.suppress_clique_check:
            final_cliques = []
            with open(os.path.join(in_args.outdir, "cliques.log"), "a") as clique_log_file:
                clique_log_file.write("""\
   #################################
   #  Initiating Clique Detection  #
   #################################

""")
                for clstr in sorted(final_clusters, key=lambda x: len(x.seq_ids), reverse=True):
                    final_cliques += clstr.create_rbh_cliques(clique_log_file)
            final_clusters = final_cliques

        # Fold singletons and doublets back into groups.
        if not in_args.suppress_singlet_folding:
            orphans = Orphans(seqbuddy=sequences, clusters=final_clusters,
                              sql_broker=broker, psi_pred_ss2_dfs=psi_pred_files)
            orphans.place_orphans()
            final_clusters = orphans.clusters
            with open(os.path.join(in_args.outdir, "orphans.log"), "a") as orphan_log_file:
                orphan_log_file.write("""\
   #################################
   #  Initiating Orphan Placement  #
   #################################

""")
                orphan_log_file.write(orphans.tmp_file.read())

        if in_args.suppress_iteration:
            break

    if not in_args.suppress_clique_check or not in_args.suppress_singlet_folding:
        logging.warning("\t-- finished in %s --" % TIMER.split())

    logging.warning("\nPlacing any collapsed paralogs into their respective clusters")
    for clust in final_clusters:
        if clust.collapsed_genes:
            with open(os.path.join(in_args.outdir, "paralog_cliques"), "a") as outfile:
                outfile.write("# %s\n" % clust.name())
                json.dump(clust.collapsed_genes, outfile)
                outfile.write("\n\n")
            for gene_id, paralog_list in clust.collapsed_genes.items():
                clust.seq_ids += paralog_list
                clust.taxa[gene_id.split(in_args.taxa_sep)[0]].append(gene_id)
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

    if in_args.workdb:
        heartbeat.end()
    broker.close()


def check_mafft_version():
    mafft_version = Popen("%s --version" % MAFFT, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    mafft_version = mafft_version[1].decode()
    mafft_version = float(re.search("v([0-9]+\.[0-9]+)", mafft_version).group(1))
    return mafft_version


def setup():
    global MAFFT
    sys.stdout.write("\033[1mWelcome to RD-MCL!\033[m\nConfirming installation...\n\n")

    sys.stdout.write("\033[1mChecking for MAFFT:\033[m ")
    try_install = False
    if os.path.isfile(MAFFT):
        ver = check_mafft_version()
        if ver < 7.245:
            sys.stdout.write("\033[91mVersion out of date (%s)\033[39m\n\n" % ver)
            try_install = True
    else:
        sys.stdout.write("\033[91mMissing\033[39m\n\n")
        try_install = True

    if try_install:
        if br.ask("Would you like the setup script to try and install MAFFT? [y]/n:"):
            if shutil.which("conda"):
                sys.stdout.write("\033[1mCalling conda...\033[m\n")
                Popen("conda install -y -c biocore mafft", shell=True).wait()
            elif shutil.which("wget") and shutil.which("make"):
                if os.path.isdir("%s%smafft" % (SCRIPT_PATH, os.sep)):
                    shutil.rmtree("%s%smafft" % (SCRIPT_PATH, os.sep))
                os.makedirs("{0}{1}mafft".format(SCRIPT_PATH, os.sep))

                cwd = os.getcwd()
                tmp_dir = br.TempDir()
                os.chdir(tmp_dir.path)
                url = "http://mafft.cbrc.jp/alignment/software/mafft-7.309-without-extensions-src.tgz"
                sys.stdout.write("\n\033[1mDownloading source code from %s\033[m\n\n" % url)
                Popen("wget %s" % url, shell=True).wait()

                sys.stdout.write("\n\033[1mUnpacking...\033[m\n")
                Popen("tar -xzf mafft-7.309-without-extensions-src.tgz", shell=True).wait()

                sys.stdout.write("\n\033[1mBuilding...\033[m\n")
                os.chdir("mafft-7.309-without-extensions" + os.sep + "core")
                with open("Makefile", "r") as ifile:
                    makefile = ifile.read()

                makefile = re.sub("PREFIX = /usr/local", "PREFIX = {0}{1}mafft".format(SCRIPT_PATH, os.sep),
                                  makefile)
                with open("Makefile", "w") as ofile:
                    ofile.write(makefile)
                Popen("make clean; make; make install", shell=True).wait()
                os.chdir(cwd)

        MAFFT = shutil.which("mafft") if shutil.which("mafft") \
            else os.path.join(SCRIPT_PATH, "mafft", "bin", "mafft")

        if not os.path.isfile(MAFFT) or check_mafft_version() < 7.245:
            sys.stdout.write("\033[91mFailed to install MAFFT.\033[39m\nPlease install the software yourself from "
                             "http://mafft.cbrc.jp/alignment/software/\n\n")
            return
    else:
        sys.stdout.write("\033[92mFound\033[39m\n")

    sys.stdout.write("\033[1mChecking for PSIPRED:\033[m ")
    path_install = []
    local_install = []
    programs = ["seq2mtx", "psipred", "psipass2"]
    for prog in programs:
        if shutil.which(prog):
            path_install.append(prog)
        elif os.path.isfile("{0}{1}psipred{1}bin{1}{2}".format(SCRIPT_PATH, os.sep, prog)):
            local_install.append(prog)
    if path_install == programs:
        sys.stdout.write("\033[92mFound\033[39m\n")
        path_install, local_install = True, False
    elif local_install == programs:
        sys.stdout.write("\033[92mFound\033[39m\n")
        path_install, local_install = False, True
    else:
        sys.stdout.write("\033[91mMissing\033[39m\n\n")
        if br.ask("Would you like the setup script to try and install PSIPRED? [y]/n:"):
            if shutil.which("conda"):
                sys.stdout.write("\033[1mCalling conda...\033[m\n")
                path_install, local_install = True, False
                Popen("conda install -y -c biocore psipred", shell=True).wait()
            elif shutil.which("wget"):
                path_install, local_install = False, True
                cwd = os.getcwd()
                tmp_dir = br.TempDir()
                os.chdir(tmp_dir.path)
                if sys.platform == "darwin":
                    version = "osx-64"
                else:
                    version = "linux-64"
                url = "https://anaconda.org/biocore/psipred/4.01/download/%s/psipred-4.01-1.tar.bz2" % version
                sys.stdout.write("\n\033[1mDownloading %s binaries from %s\033[m\n\n" % (version, url))
                Popen("wget %s" % url, shell=True).wait()

                sys.stdout.write("\n\033[1mUnpacking...\033[m\n")
                Popen("tar -xjf psipred-4.01-1.tar.bz2", shell=True).wait()

                sys.stdout.write("\n\033[1mInstalling...\033[m\n")
                if os.path.isdir("%s%spsipred" % (SCRIPT_PATH, os.sep)):
                    shutil.rmtree("%s%spsipred" % (SCRIPT_PATH, os.sep))
                os.makedirs("%s%spsipred" % (SCRIPT_PATH, os.sep))

                shutil.move("bin", "{0}{1}psipred{1}".format(SCRIPT_PATH, os.sep))
                shutil.move("share{0}psipred_4.01{0}data".format(os.sep),
                            "{0}{1}psipred{1}".format(SCRIPT_PATH, os.sep))
                os.chdir(cwd)
            else:
                sys.stdout.write("\033[91mFailed to install PSIPRED.\033[39m\nPlease install the software yourself from"
                                 " http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/\n\n")
                return
    # Confirm all psipred weight files are in the rdmcl directory
    weight_files = ["weights.dat", "weights.dat2", "weights.dat3", "weights_p2.dat",
                    "weights_s.dat", "weights_s.dat2", "weights_s.dat3"]
    error_msg = """\033[1mError:\033[m psi-pred data file '{0}' not found in {1}!
    Please try reinstalling PSIPRED:

       $: conda install -c biocore --no-deps --force psipred

    or build from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/

    If the problem persists, please create an issue at https://github.com/biologyguy/RD-MCL/issues
    """

    data_dir = "{0}{1}psipred{1}data".format(SCRIPT_PATH, os.sep)

    if path_install:
        os.makedirs("%s%spsipred" % (SCRIPT_PATH, os.sep), exist_ok=True)
        os.makedirs(data_dir, exist_ok=True)

        psipred_bin_dir = shutil.which("psipred").split(os.sep)[:-2]
        root, dirs, files = next(os.walk(os.sep + os.path.join(*psipred_bin_dir, "share")))
        psipred_data_dir = re.search("'(psipred.*?)'[,\]]", str(dirs)).group(1)
        psipred_data_dir = os.sep + os.path.join(*psipred_bin_dir, "share", psipred_data_dir, "data")
        for next_file in weight_files:
            if not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)) \
                    and not os.path.isfile(os.path.join(psipred_data_dir, next_file)):
                print(error_msg.format(next_file, psipred_data_dir))
                return
            elif not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)):
                shutil.copyfile(os.path.join(psipred_data_dir, next_file),
                                "{0}{1}{2}".format(data_dir, os.sep, next_file))
    elif local_install:
        for next_file in weight_files:
            if not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)):
                print(error_msg.format(next_file, data_dir))
                return

    open(os.path.join(SCRIPT_PATH, "config.ini"), "w").close()
    sys.stdout.write("\n\033[1mSuccess! You're all set.\033[m\n")
    return


def main():
    in_args = argparse_init()
    if 'setup' in in_args:
        setup()
    else:
        full_run(in_args)
    return


if __name__ == '__main__':
    '''
    cteno_panxs = Sb.SeqBuddy("tests/unit_test_resources/Cteno_pannexins.fa")
    ids = [rec.id for rec in cteno_panxs.records]
    sim_scores = pd.read_csv("tests/unit_test_resources/Cteno_pannexins_sim.scores", "\t", index_col=False, header=None)
    sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]

    cluster = Cluster(ids, sim_scores)
    cluster.score()
    '''
    main()
