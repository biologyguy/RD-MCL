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
from os.path import join
import re
import shutil
import json
import logging
import time
import argparse
import sqlite3
from io import StringIO
from subprocess import Popen, PIPE
from multiprocessing import Lock, Process
from random import choice, Random, randint, random
from math import log2
from collections import OrderedDict
from copy import deepcopy, copy
# from hashlib import md5


# 3rd party
import pandas as pd
import numpy as np
# import statsmodels.api as sm
# import statsmodels.stats.api as sms
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br

# Other RD-MCL modules
try:
    import mcmcmc
    import helpers as hlp
    from install import setup
except ImportError:
    from . import mcmcmc
    from . import helpers as hlp
    from .install import setup


# Globals (I apologize to all of you cringing over this...)
SCRIPT_PATH = hlp.SCRIPT_PATH
VERSION = hlp.VERSION
VERSION.name = "rdmcl"
ALIGNMETHOD = "clustalo"
ALIGNPARAMS = ""
LOCK = Lock()
MULTICORE_LOCK = Lock()
PROGRESS_LOCK = Lock()
CPUS = br.usable_cpu_count()
TIMER = hlp.Timer()
MIN_SIZE_TO_WORKER = 15  # (MIN_SIZE_TO_WORKER**2 - MIN_SIZE_TO_WORKER) / 2  =  105
GAP_OPEN = -5
GAP_EXTEND = 0
BLOSUM62 = hlp.make_full_mat(SeqMat(MatrixInfo.blosum62))
MCMC_CHAINS = 3
GELMAN_RUBIN = 1.1
WORKER_OUT = ""
WORKER_DB = ""
HEARTBEAT_DB = ""
MAX_WORKER_WAIT = 240
MASTER_ID = None
MASTER_PULSE = 60
PSIPREDDIR = ""
TRIMAL = ["gappyout", 0.5, 0.75, 0.9, 0.95, "clean"]

if os.path.isfile(join(SCRIPT_PATH, "hmmer", "hmm_fwd_back")):
    HMM_FWD_BACK = join(SCRIPT_PATH, "hmmer", "hmm_fwd_back")
elif shutil.which("hmm_fwd_back"):
    HMM_FWD_BACK = shutil.which("hmm_fwd_back")
elif "-setup" in sys.argv:
    HMM_FWD_BACK = ""
else:
    sys.stderr.write("Error: hmm_fwd_back program not found. Please run `rdmcl.py -setup` to fix this.\n")
    sys.exit()

if os.path.isfile(join(SCRIPT_PATH, "hmmer", "hmmbuild")):
    HMMBUILD = join(SCRIPT_PATH, "hmmer", "hmmbuild")
elif shutil.which("hmmbuild"):
    HMMBUILD = shutil.which("hmmbuild")
elif "-setup" in sys.argv:
    HMM_FWD_BACK = ""
else:
    sys.stderr.write("Error: hmmbuild program not found. Please run `rdmcl.py -setup` to fix this.\n")
    sys.exit()

ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
               "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
for aa in ambiguous_X:
    pair = sorted((aa, "X"))
    pair = tuple(pair)
    BLOSUM62[pair] = ambiguous_X[aa]

# Set global precision levels to prevent weird rounding on different architectures
np.set_printoptions(precision=12)
pd.set_option("display.precision", 12)


# ToDo: Maybe remove support for taxa_sep. It's complicating my life, so just impose the '-' on users?
class Cluster(object):
    def __init__(self, seq_ids, sim_scores, taxa_sep="-", group_prefix="group",
                 parent=None, collapse=False, r_seed=None, fwd_scores=None):
        """
        - Note that reciprocal best hits between paralogs are collapsed when instantiating group_0, so
          no problem strongly penalizing all paralogs in the scoring algorithm

        :param seq_ids: Sequence IDs
        :type seq_ids: list or set
        :param sim_scores: All-by-all similarity matrix for everything in the sequence_ids (and parental clusters)
        :type sim_scores: pandas.DataFrame
        :param taxa_sep: What character partitions taxon from rec id (deprecated)
        :param group_prefix: The name that is attached to each group
        :param parent: Parental sequence_ids
        :type parent: Cluster
        :param collapse: Specify whether reciprocal best hit paralogs should be combined
        :type collapse: bool
        :param r_seed: Set the random generator seed value
        """
        # Do an initial sanity check on the incoming graph and list of sequence ids
        ids_len = len(seq_ids)
        if len(set(seq_ids)) != ids_len:
            raise AttributeError("seq_ids are not all unique.")
        seq_ids = sorted(seq_ids)

        expected_num_edges = int(((ids_len**2) - ids_len) / 2)
        if len(sim_scores.index) != expected_num_edges:
            seq_id_hash = hlp.md5_hash(str(", ".join(seq_ids)))
            sim_scores_ids = sim_scores.seq1.append(sim_scores.seq2)
            sim_scores_ids = sorted(list(sim_scores_ids.unique()))

            with open("%s.error" % seq_id_hash, "w") as ofile:
                ofile.write("Parent: %s\nCollapse: %s\n\nSeq IDs\n%s\n\nSimScore IDs\n%s" %
                            (parent, collapse, seq_ids, sim_scores_ids))
            raise ValueError("The number of incoming sequence ids (%s) does not match the expected graph size of %s"
                             " rows (observed %s rows).\nError report at %s.error" %
                             (ids_len, expected_num_edges, len(sim_scores.index), seq_id_hash))

        self.sim_scores = sim_scores
        self.taxa_sep = taxa_sep
        self.parent = parent

        self.subgroup_counter = 0
        self.cluster_score = None
        self.collapsed_genes = OrderedDict()  # If paralogs are reciprocal best hits, collapse them
        self.rand_gen = Random(r_seed)
        self.fwd_scores = fwd_scores
        self._name = None
        self.taxa = OrderedDict()  # key = taxa id. value = list of genes coming fom that taxa.

        for next_seq_id in seq_ids:
            taxa = next_seq_id.split(taxa_sep)[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)

        self.seq_ids = set(seq_ids)
        if parent:
            self.max_genes_in_a_taxa = parent.max_genes_in_a_taxa
            for representative, collapsed_genes in parent.collapsed_genes.items():
                if representative in seq_ids:
                    self.collapsed_genes[representative] = collapsed_genes

        else:
            self.max_genes_in_a_taxa = max([len(self.taxa[taxa]) for taxa in self.taxa])
            self._name = "%s_0" % group_prefix
            # This next bit can collapse all paralog reciprocal best-hit cliques so they don't gum up MCL
            if collapse:
                column = "r_square" if "r_square" in self.sim_scores.columns else "raw_score"
                self.collapse(column)

        self.seq_ids_str = str(", ".join(sorted(self.seq_ids)))
        self.seq_id_hash = hlp.md5_hash(self.seq_ids_str)

    def reset_seq_ids(self, seq_ids):
        # Note that this does NOT reset sim_scores. This needs to be updated manually
        seq_ids = sorted(seq_ids)
        self.seq_ids_str = str(", ".join(seq_ids))
        self.seq_id_hash = hlp.md5_hash(self.seq_ids_str)
        self.seq_ids = set(seq_ids)
        self.taxa = OrderedDict()
        for next_seq_id in self.seq_ids:
            taxa = next_seq_id.split(self.taxa_sep)[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)

    def collapse(self, column):
        breakout = False
        seq_ids = sorted(self.seq_ids)
        while not breakout:
            breakout = True
            for seq1_id in seq_ids:
                seq1_taxa = seq1_id.split(self.taxa_sep)[0]
                paralog_best_hits = []
                for best_hits_seq1 in self.get_best_hits(seq1_id, column).itertuples():  # grab best hit for seq1
                    seq2_id = best_hits_seq1.seq1 if best_hits_seq1.seq1 != seq1_id else best_hits_seq1.seq2
                    if seq2_id.split(self.taxa_sep)[0] != seq1_taxa:
                        paralog_best_hits = []
                        break
                    for best_hits_seq2 in self.get_best_hits(seq2_id, column).itertuples():  # Confirm reciprocal
                        if seq1_id in [best_hits_seq2.seq1, best_hits_seq2.seq2]:
                            paralog_best_hits.append(seq2_id)
                            break

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
        self.seq_ids = set(seq_ids)
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

    def get_best_hits(self, gene, column):
        """
        Search through the sim scores data frame to find the gene(s) with the highest similarity to the query
        :param gene: Query gene name
        :type gene: str
        :param column: Which column in self.sim_scores are we assessing?
        :return: The best hit(s). Multiple matches will all have the same similarity score.
        :rtype: pd.DataFrame
        """
        best_hits = self.sim_scores[(self.sim_scores.seq1 == gene) | (self.sim_scores.seq2 == gene)]
        if not best_hits.empty:
            best_hits = best_hits.loc[best_hits[column] == max(best_hits[column])].values
            best_hits = pd.DataFrame(best_hits, columns=list(self.sim_scores.columns.values))
        return best_hits

    def recursive_best_hits(self, gene, global_best_hits, tested_ids, column):
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
        :param column: Which column in self.sim_scores are we assessing?
        :return: The final clique
        :rtype: pd.DataFrame
        """
        best_hits = self.get_best_hits(gene, column)
        global_best_hits = global_best_hits.append(best_hits, ignore_index=True)
        for _edge in best_hits.itertuples():
            if _edge.seq1 not in tested_ids:
                tested_ids.append(_edge.seq1)
                global_best_hits = self.recursive_best_hits(_edge.seq1, global_best_hits, tested_ids, column)
            if _edge.seq2 not in tested_ids:
                tested_ids.append(_edge.seq2)
                global_best_hits = self.recursive_best_hits(_edge.seq2, global_best_hits, tested_ids, column)
        return global_best_hits.drop_duplicates()

    def perturb(self, scores, col_name="score"):
        """
        Add a very small bit of variance into score values when std == 0
        :param scores: Similarity scores
        :type scores: pd.DataFrame
        :param col_name: DataFrame column label
        :return: updated data frame
        """
        valve = br.SafetyValve(global_reps=10)
        while scores[col_name].std() == 0:
            valve.step("Failed to perturb:\n%s" % scores)
            for indx, score in scores[col_name].iteritems():
                scores.at[indx, col_name] = round(self.rand_gen.gauss(score, (score * 0.0000001)), 12)
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
            base_cluster.seq_ids.add(orphan)
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

        dim_ret_base_score = self.get_dim_ret_base_score()
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
            subscore *= dim_ret_base_score ** indx
            score += subscore

        self.cluster_score = score
        return self.cluster_score

    def get_dim_ret_base_score(self):
        ave_num_paralogs = len(self.seq_ids) / len(self.taxa)
        if ave_num_paralogs < len(self.taxa):  # If more taxa than paralogs => DRB < 0.5, easier
            return (ave_num_paralogs / len(self.taxa)) * 0.5
        elif ave_num_paralogs == len(self.taxa):  # If num taxa == num paralogs => DRB = 0.5
            return 0.5
        else:  # If more paralogs than taxa => DRB > 0.5
            return 1 - ((len(self.taxa) / ave_num_paralogs) * 0.5)

    def _score_direct_replicate_penalty(self):
        # Final scores can be negative, which is problematic
        # This is currently only accessible by modifying the default in score()
        unique_taxa = 0
        replicate_taxa = 0
        base_cluster = deepcopy(self.parent) if self.parent else self

        # It's possible that a cluster can have sequences not present in the parent (i.e., following orphan placement);
        # check for this, and add the sequences to base_cluster if necessary.
        items_not_in_parent = set(self.seq_ids) - set(base_cluster.seq_ids)
        for orphan in items_not_in_parent:
            base_cluster.seq_ids.add(orphan)
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
        log_file.write("# ####### Testing %s ####### #\n%s\n" % (self.name(), sorted(self.seq_ids)))
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
                                                    [gene], "raw_score")
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
                                             & (outer_scores.seq1.isin(clique_ids) == False))]   # Note the == must stay

                # if all sim scores in a group are identical, we can't get a KDE. Fix by perturbing the scores a little.
                clique_scores = self.perturb(clique_scores)
                clique_scores = self.perturb(clique_scores, "raw_score")

                outer_scores = self.perturb(outer_scores)
                outer_scores = self.perturb(outer_scores, "raw_score")

                total_kde = scipy.stats.gaussian_kde(outer_scores.raw_score, bw_method='silverman')
                log_file.write("\t\t\tOuter KDE: {'shape': %s, 'covariance': %s, 'inv_cov': %s, '_norm_factor': %s}\n" %
                               (total_kde.dataset.shape, round(total_kde.covariance[0][0], 12),
                                round(total_kde.inv_cov[0][0], 12), round(total_kde._norm_factor, 12)))

                clique_kde = scipy.stats.gaussian_kde(clique_scores.raw_score, bw_method='silverman')
                log_file.write("\t\t\tClique KDE: {'shape': %s, 'covariance': %s, "
                               "'inv_cov': %s,  '_norm_factor': %s}\n" %
                               (clique_kde.dataset.shape, round(clique_kde.covariance[0][0], 12),
                                round(clique_kde.inv_cov[0][0], 12), round(clique_kde._norm_factor, 12)))

                clique_resample = clique_kde.resample(10000)[0]  # ToDo: figure out how to control this with r_seed!
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
        remaining_seqs = sorted(self.seq_ids)
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
                       "\n\t\t".join([str(sorted(res.seq_ids)) for res in results]))
        return results

    def __len__(self):
        return len(self.seq_ids)

    def __str__(self):
        return str(sorted(self.seq_ids))


def cluster2database(cluster, sql_broker, alignment):
    """
    Update the database with a cluster
    :param cluster: Cluster object
    :param sql_broker: A running helpers.SQLiteBroker object
    :param alignment: An alignment as an AlignBuddy object or string
    :return:
    """
    sql_broker.query("""INSERT OR IGNORE INTO data_table (hash, seq_ids, alignment, graph, cluster_score)
                        VALUES (?, ?, ?, ?, ?)
                        """, (cluster.seq_id_hash, cluster.seq_ids_str, str(alignment),
                              cluster.sim_scores.to_csv(header=None, index=False), cluster.score(),))
    return


# ################ PSI-PRED FUNCTIONS ################ #
def mc_psi_pred(seq_obj, args):
    outdir = args[0]
    if os.path.isfile(join(outdir, "%s.ss2" % seq_obj.id)):
        return
    result = run_psi_pred(seq_obj)
    with open(join(outdir, "%s.ss2" % seq_obj.id), "w") as ofile:
        ofile.write(result)
    return


def run_psi_pred(seq_rec):
    temp_dir = br.TempDir()
    pwd = os.getcwd()
    psipred_dir = join(SCRIPT_PATH, "psipred")
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
        data_weights = join(psipred_dir, "data", "weights")
        command = '''\
    {0}{3}bin{3}seq2mtx sequence.fa > {1}{3}{2}.mtx;
    {0}{3}bin{3}psipred {1}{3}{2}.mtx {4}.dat {4}.dat2 {4}.dat3 > {1}{3}{2}.ss;
    {0}{3}bin{3}psipass2 {4}_p2.dat 1 1.0 1.0 {1}{3}{2}.ss2 {1}{3}{2}.ss > {1}{3}{2}.horiz;
    '''.format(psipred_dir, temp_dir.path, seq_rec.id, os.sep, data_weights)

    Popen(command, shell=True).wait()
    os.chdir(pwd)
    with open(join(temp_dir.path, "%s.ss2" % seq_rec.id), "r") as ifile:
        result = ifile.read()
    return result


def read_ss2_file(path):
    ss_file = pd.read_csv(path, comment="#", header=None, delim_whitespace=True)
    ss_file.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
    return ss_file


def compare_psi_pred(psi1_df, psi2_df):
    num_extra_gaps = 0
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
            num_extra_gaps += 1
    align_len = len(psi2_df) + num_extra_gaps  # Note that any gaps in psi1 are automatically accounted for by len(psi2)
    ss_score /= align_len
    return ss_score
# ################ END PSI-PRED FUNCTIONS ################ #


def orthogroup_caller(master_cluster, cluster_list, seqbuddy, sql_broker, progress, outdir, psi_pred_ss2,
                      steps=1000, chains=3, walkers=2, quiet=True, taxa_sep="-", r_seed=None, convergence=None,
                      resume=False):
    """
    Run MCMCMC on MCL to find the best orthogroups
    :param master_cluster: The group to be subdivided
    :type master_cluster: Cluster
    :param cluster_list: When a sequence_ids is finalized after recursion, it is appended to this list
    :param seqbuddy: The sequences that are included in the master sequence_ids
    :param sql_broker: Multithread SQL broker that can be queried
    :param progress: Progress class
    :param outdir: where are files being written to?
    :param psi_pred_ss2: OrdredDict of all ss2 dataframes paths with record IDs as key
    :param steps: How many MCMCMC iterations to run TODO: calculate this on the fly
    :param chains: Number of MCMCMC chains to spin off
    :param walkers: Number of Metropolis-Hastings walkers per chain
    :param quiet: Suppress StdErr
    :param taxa_sep: The string that separates taxon names from gene names
    :param r_seed: Set the random generator seed value
    :param convergence: Set minimum Gelman-Rubin PSRF value for convergence
    :param resume: Try to pick up from a previous run
    :return: list of sequence_ids objects
    """
    def save_cluster(end_message=None):
        cluster_list.append(master_cluster)
        if not os.path.isfile(join(mcmcmc_path, "best_group")):
            with open(join(mcmcmc_path, "best_group"), "w") as _ofile:
                _ofile.write('\t'.join(master_cluster.seq_ids))
        if end_message:
            with open(join(mcmcmc_path, "end_message.log"), "w") as _ofile:
                _ofile.write(end_message + "\n")

        _, alignment = retrieve_all_by_all_scores(seqbuddy, psi_pred_ss2, sql_broker, quiet=True)
        alignment.write(join(outdir, "alignments", master_cluster.name()))
        master_cluster.sim_scores.to_csv(join(outdir, "sim_scores", "%s.scores" % master_cluster.name()),
                                         header=None, index=False, sep="\t")
        alignment = Alb.generate_hmm(alignment, HMMBUILD)
        with open(join(outdir, "hmm", master_cluster.name()), "w") as _ofile:
            _ofile.write(alignment.alignments[0].hmm)
        update = len(master_cluster.seq_ids) if not master_cluster.subgroup_counter else 0
        progress.update("placed", update)
        return

    rand_gen = Random(r_seed)
    master_cluster.set_name()
    mcmcmc_path = join(outdir, "mcmcmc", master_cluster.name())
    os.makedirs(mcmcmc_path, exist_ok=True)
    open(join(mcmcmc_path, "max.txt"), "w").close()
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

    # I know what the best and worst possible scores are, so let MCMCMC know (better for calculating acceptance rates)
    # The worst score possible would be all genes in each taxa segregated.
    worst_possible_score = 0
    for taxon, seq_ids in master_cluster.taxa.items():
        bad_clust = Cluster(seq_ids, master_cluster.pull_scores_subgraph(seq_ids), parent=master_cluster)
        worst_possible_score += bad_clust.score()

    # The best score would be perfect separation of taxa
    sub_clusters = [[]]
    for taxon, genes in master_cluster.taxa.items():
        for indx, gene in enumerate(genes):
            if len(sub_clusters) > indx:
                sub_clusters[indx].append(gene)
            else:
                sub_clusters.append([gene])
    best_possible_score = 0
    for seq_ids in sub_clusters:
        subcluster = Cluster(seq_ids, master_cluster.pull_scores_subgraph(seq_ids), parent=master_cluster)
        best_possible_score += subcluster.score()

    if best_possible_score == worst_possible_score:
        return cluster_list

    mcmcmc_params = [mcmcmc_path, seqbuddy, master_cluster,
                     taxa_sep, sql_broker, psi_pred_ss2, progress, chains * (walkers + 2)]
    mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1, quiet=quiet,
                                   num_walkers=walkers, num_chains=chains, convergence=convergence,
                                   outfile_root=join(mcmcmc_path, "mcmcmc_out"), params=mcmcmc_params,
                                   include_lava=True, include_ice=True, r_seed=rand_gen.randint(1, 999999999999999),
                                   min_max=(worst_possible_score, best_possible_score))

    mcmcmc_factory.reset_params([mcmcmc_path, seqbuddy, master_cluster,
                                 taxa_sep, sql_broker, psi_pred_ss2, progress, chains * (walkers + 2)])

    if resume:
        if not mcmcmc_factory.resume():
            mcmcmc_factory.run()

    else:
        mcmcmc_factory.run()

    best_score = pd.DataFrame()
    for indx in range(len(mcmcmc_factory.chains)):
        mcmcmc_output = pd.read_csv(join(mcmcmc_path, "mcmcmc_out_%s.csv" % (indx+1)), "\t", index_col=False)
        result = mcmcmc_output.loc[mcmcmc_output["result"] == mcmcmc_output["result"].max()]
        best_score = result if best_score.empty or result["result"].iloc[0] > best_score["result"].iloc[0] \
            else best_score

    if round(best_score["result"].iloc[0], 8) <= round(master_cluster.score(), 8):
        save_cluster("New best score of %s is ≤ master cluster at %s"
                     % (round(best_score["result"].iloc[0], 8), round(master_cluster.score(), 8)))
        return cluster_list

    mcl_obj = hlp.MarkovClustering(master_cluster.sim_scores, inflation=best_score["I"].iloc[0],
                                   edge_sim_threshold=best_score["gq"].iloc[0])
    mcl_obj.run()
    progress.update('mcl_runs', 1)
    mcl_clusters = mcl_obj.clusters

    # Write out the actual best clusters
    best_clusters = ['\t'.join(cluster) for cluster in mcl_clusters]
    with open(join(mcmcmc_path, "best_group"), "w") as ofile:
        ofile.write('\n'.join(best_clusters))

    recursion_clusters = []
    for sub_cluster in mcl_clusters:
        cluster_ids_hash = hlp.md5_hash(", ".join(sorted(sub_cluster)))
        if len(sub_cluster) == 1:
            sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
        else:
            # All mcl sub clusters are written to database in mcmcmc_mcl(), so no need to check if exists
            graph = sql_broker.query("SELECT (graph) FROM data_table WHERE hash=?", (cluster_ids_hash,))[0][0]

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
            _, align = retrieve_all_by_all_scores(Sb.pull_recs(Sb.make_copy(seqbuddy),
                                                               "^%s$" % ("$|^".join(sub_cluster.seq_ids))),
                                                  psi_pred_ss2, sql_broker, quiet=True)
            align = Alb.generate_hmm(align, HMMBUILD)
            with open(join(outdir, "hmm", sub_cluster.name()), "w") as ofile:
                ofile.write(align.alignments[0].hmm)
            cluster_list.append(sub_cluster)
            continue
        recursion_clusters.append(sub_cluster)

    for sub_cluster in recursion_clusters:
        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in sub_cluster.seq_ids])

        # Recursion... Reassign cluster_list, as all clusters are returned at the end of a call to orthogroup_caller
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, sql_broker=sql_broker,
                                         progress=progress, outdir=outdir, steps=steps, quiet=quiet, chains=chains,
                                         walkers=walkers, taxa_sep=taxa_sep, convergence=convergence, resume=resume,
                                         r_seed=rand_gen.randint(1, 999999999999999), psi_pred_ss2=psi_pred_ss2)

    save_cluster("Sub clusters returned")
    return cluster_list


# #########  Miscellaneous  ########## #
class Progress(object):
    def __init__(self, outdir, base_cluster):
        self.outdir = outdir
        with open(join(self.outdir, ".progress"), "w") as progress_file:
            _progress = {"mcl_runs": 0, "placed": 0, "total": len(base_cluster)}
            json.dump(_progress, progress_file)

    def update(self, key, value):
        with PROGRESS_LOCK:
            with open(join(self.outdir, ".progress"), "r") as ifile:
                _progress = json.load(ifile)
                _progress[key] += value
            with open(join(self.outdir, ".progress"), "w") as _ofile:
                json.dump(_progress, _ofile)
        return

    def read(self):
        with PROGRESS_LOCK:
            with open(join(self.outdir, ".progress"), "r") as ifile:
                return json.load(ifile)

    def __str__(self):
        _progress = self.read()
        return "MCL runs processed: %s. Sequences placed: %s/%s. Run time: " \
               % (_progress['mcl_runs'], _progress['placed'], _progress['total'])


def check_sequences(seqbuddy, taxa_sep):
    logging.warning("Checking that the format of all sequence ids matches 'taxa%sgene'" % taxa_sep)
    failures = []
    taxa = []
    for rec in seqbuddy.records:
        rec_id = rec.id.split(taxa_sep)
        if len(rec_id) != 2:
            failures.append(rec.id)
        else:
            taxa.append(rec_id[0])
    if failures:
        logging.error("Malformed sequence id(s): '%s'\nThe taxa separator character is currently set to '%s',\n"
                      " which can be changed with the '-ts' flag" % (", ".join(failures), taxa_sep))
        return False
    else:
        logging.warning("    %s sequences PASSED" % len(seqbuddy))
        logging.info("%s unique taxa present" % len(set(taxa)))
        return True


class HeartBeat(object):
    def __init__(self, hbdb_path, pulse_rate, thread_type='master', dummy=False):
        self.hbdb_path = hbdb_path
        self.pulse_rate = pulse_rate
        self.id = None
        self.running_process = None
        self.check_file = br.TempFile()
        self.thread_type = thread_type
        self.dummy = dummy  # This allows an object with start() and end() methods that won't spin off a daemon
        if not dummy:
            with hlp.ExclusiveConnect(self.hbdb_path) as cursor:
                try:
                    cursor.execute('CREATE TABLE heartbeat (thread_id INTEGER PRIMARY KEY AUTOINCREMENT, '
                                   'thread_type TEXT, pulse INTEGER)')
                except sqlite3.OperationalError:
                    pass

    def _run(self, check_file_path):
        if self.dummy:
            return
        split_time = time.time()
        try:
            while True:
                with open(check_file_path, "r") as ifile:
                    ifile_content = ifile.read()
                if ifile_content != "Running":
                    break

                if split_time < time.time() - self.pulse_rate:
                    with hlp.ExclusiveConnect(self.hbdb_path) as cursor:
                        num_threads = cursor.execute("SELECT COUNT(*) FROM heartbeat").fetchone()[0]
                        activity_modifier = num_threads * 0.1
                        cursor.execute("INSERT or REPLACE INTO heartbeat (thread_id, thread_type, pulse) VALUES"
                                       " (?, ?, ?)", (self.id, self.thread_type,
                                                      round(time.time() + cursor.lag + activity_modifier),))
                    split_time = time.time() + activity_modifier
                time.sleep(random())
        except KeyboardInterrupt:
            open(check_file_path, "w").close()
        return

    def start(self):
        if self.dummy:
            return
        if self.running_process:
            self.end()
        self.check_file.write("Running")
        with hlp.ExclusiveConnect(self.hbdb_path) as cursor:
            cursor.execute("INSERT INTO heartbeat (thread_type, pulse) "
                           "VALUES (?, ?)", (self.thread_type, round(time.time()),))
            self.id = cursor.lastrowid
        p = Process(target=self._run, args=(self.check_file.path,))
        p.daemon = 1
        p.start()
        self.running_process = p
        return

    def end(self):
        if self.dummy:
            return
        if not self.running_process:
            return
        self.check_file.clear()
        while self.running_process.is_alive():
            continue
        self.running_process = None
        with hlp.ExclusiveConnect(self.hbdb_path) as cursor:
            cursor.execute("DELETE FROM heartbeat WHERE thread_id=?", (self.id,))
        self.id = None
        return


# ################ SCORING FUNCTIONS ################ #
def mc_score_sequences(seq_pairs, args):
    # ##################################################################### #
    # Calculate the best possible scores, and divide by the observed scores #
    # ##################################################################### #

    results = ["" for _ in seq_pairs]
    alb_obj, gap_open, gap_extend, output_file = args
    for indx, seq_pair in enumerate(seq_pairs):
        id1, id2, ss2df1, ss2df2 = seq_pair
        id_regex = "^%s$|^%s$" % (id1, id2)
        # Alignment comparison
        alb_copy = Alb.make_copy(alb_obj)
        Alb.pull_records(alb_copy, id_regex)

        subs_mat_score = compare_pairwise_alignment(alb_copy, gap_open, gap_extend)

        # PSI PRED comparison
        ss_score = compare_psi_pred(ss2df1, ss2df2)
        results[indx] = "\n%s,%s,%s,%s" % (id1, id2, subs_mat_score, ss_score)
    with LOCK:
        with open(output_file, "a") as ofile:
            ofile.write("".join(results))
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


def mc_create_all_by_all_scores(seqbuddy, args):
    psi_pred_ss2, sql_broker = args
    retrieve_all_by_all_scores(seqbuddy, psi_pred_ss2, sql_broker, quiet=True)
    return


def retrieve_all_by_all_scores(seqbuddy, psi_pred_ss2, sql_broker, quiet=False):
    """
    :param seqbuddy: SeqBuddy object
    :param psi_pred_ss2: OrderedDict of {seqID: ss2 dataframe path}
    :param sql_broker: Active broker object to search/update SQL database
    :param quiet: Supress multicore output
    :return: sim_scores, Alb.AlignBuddy
    """
    seq_ids = sorted([rec.id for rec in seqbuddy.records])
    seq_id_hash = hlp.md5_hash(", ".join(seq_ids))

    if len(seqbuddy) == 1:
        sim_scores = pd.DataFrame(data=None, columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
        alignment = Alb.AlignBuddy(str(seqbuddy), in_format="fasta")
        cluster2database(Cluster(seq_ids, sim_scores), sql_broker, alignment)
        return sim_scores, alignment

    # Grab from the database first, if the data exists there already
    query = sql_broker.query("SELECT graph, alignment FROM data_table WHERE hash=?", (seq_id_hash,))
    if query and len(query[0]) == 2:
        sim_scores, alignment = query[0]
        sim_scores = pd.read_csv(StringIO(sim_scores), index_col=False, header=None)
        sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
        return sim_scores, Alb.AlignBuddy(alignment, in_format="fasta")

    # Try to feed the job to independent workers
    if WORKER_DB and os.path.isfile(WORKER_DB) and len(seq_ids) > MIN_SIZE_TO_WORKER:
        workerjob = WorkerJob(seqbuddy, sql_broker)
        worker_result = workerjob.run()
        if worker_result:
            return worker_result

    # If the job is small or couldn't be pushed off on a worker, do it directly
    all_by_all_obj = AllByAllScores(seqbuddy, psi_pred_ss2, sql_broker, quiet=quiet)
    return all_by_all_obj.create()


class AllByAllScores(object):
    def __init__(self, seqbuddy, psi_pred_ss2, sql_broker, quiet=False):
        self.seqbuddy = Sb.make_copy(seqbuddy)
        self.seq_ids = sorted([rec.id for rec in self.seqbuddy.records])
        self.psi_pred_ss2 = psi_pred_ss2
        self.sql_broker = sql_broker
        self.quiet = quiet

    def create(self):
        """
        Generate a multiple sequence alignment and pull out all-by-all similarity graph
        :return:
        """
        psi_pred_ss2_dfs = OrderedDict()

        for rec in self.seqbuddy.records:
            psi_pred_ss2_dfs[rec.id] = read_ss2_file(self.psi_pred_ss2[rec.id])

        alignment = Alb.generate_msa(Sb.make_copy(self.seqbuddy), ALIGNMETHOD, ALIGNPARAMS, quiet=True)

        # Need to specify what columns the PsiPred files map to now that there are gaps.
        psi_pred_ss2_dfs = update_psipred(alignment, psi_pred_ss2_dfs, "msa")

        # Scores seem to be improved by removing gaps. ToDo: Need to test this explicitly for the paper
        # Only remove columns up to a 50% reduction in average seq length and only if all sequences are retained
        alignment = trimal(self.seqbuddy, TRIMAL, alignment)

        # Re-update PsiPred files now that some columns, possibly including non-gap characters, are removed
        psi_pred_ss2_dfs = update_psipred(alignment, psi_pred_ss2_dfs, "trimal")

        all_by_all_len, all_by_all = prepare_all_by_all(self.seqbuddy, psi_pred_ss2_dfs, CPUS)
        all_by_all_outfile = br.TempFile()
        all_by_all_outfile.write("seq1,seq2,subsmat,psi")
        score_sequences_params = [alignment, GAP_OPEN, GAP_EXTEND, all_by_all_outfile.path]
        with MULTICORE_LOCK:
            br.run_multicore_function(all_by_all, mc_score_sequences, score_sequences_params,
                                      quiet=self.quiet, max_processes=CPUS)
        sim_scores = pd.read_csv(all_by_all_outfile.get_handle("r"), index_col=False)
        sim_scores = set_final_sim_scores(sim_scores)
        cluster2database(Cluster(self.seq_ids, sim_scores), self.sql_broker, alignment)
        return sim_scores, alignment


def update_psipred(alignment, psipred_dfs, mode):
    if mode == "msa":
        for rec in alignment.records_iter():
            ss_file = psipred_dfs[rec.id]
            ss_counter = 0
            for indx, residue in enumerate(rec.seq):
                if residue != "-":
                    psipred_dfs[rec.id].at[ss_counter, "indx"] = indx
                    # psipred_dfs[rec.id].set_value(ss_counter, "indx", indx)
                    ss_counter += 1
            psipred_dfs[rec.id] = ss_file

    elif mode == "trimal":
        for rec in alignment.records_iter():
            # Instantiate list of max possible size
            new_psi_pred = [0 for _ in range(len(psipred_dfs[rec.id].index))]
            indx = 0
            for row in psipred_dfs[rec.id].itertuples():
                if alignment.alignments[0].position_map[int(row[1])][1]:
                    new_psi_pred[indx] = list(row)[1:]
                    indx += 1
            new_psi_pred = new_psi_pred[:indx]
            psipred_dfs[rec.id] = pd.DataFrame(new_psi_pred, columns=["indx", "aa", "ss", "coil_prob",
                                                                      "helix_prob", "sheet_prob"])
    else:
        raise ValueError("Unrecognized mode '%s': select from ['msa', 'trimal']" % mode)
    return psipred_dfs


def trimal(seqbuddy, trimal_modes, alignment):
        # Scores seem to be improved by removing gaps. ToDo: Need to test this explicitly for the paper
        # Only remove columns up to a 50% reduction in average seq length and only if all sequences are retained
        ave_seq_length = Sb.ave_seq_length(seqbuddy)
        for threshold in trimal_modes:
            align_copy = Alb.trimal(Alb.make_copy(alignment), threshold=threshold)
            cleaned_seqs = Sb.clean_seq(Sb.SeqBuddy(str(align_copy)))
            cleaned_seqs = Sb.delete_small(cleaned_seqs, 1)
            # Structured this way for unit test purposes
            if len(alignment.records()) != len(cleaned_seqs):
                continue
            elif Sb.ave_seq_length(cleaned_seqs) / ave_seq_length < 0.5:
                continue
            else:
                alignment = align_copy
                break
        return alignment


def prepare_all_by_all(seqbuddy, psipred_dfs, cpus):
        ids1 = [rec.id for rec in seqbuddy.records]
        ids2 = copy(ids1)
        data = [0 for _ in range(int((len(ids1)**2 - len(ids1)) / 2))]
        indx = 0
        for rec1 in ids1:
            del ids2[ids2.index(rec1)]
            for rec2 in ids2:
                data[indx] = (rec1, rec2, psipred_dfs[rec1], psipred_dfs[rec2])
                indx += 1

        return len(data), hlp.chunk_list(data, cpus)


def set_final_sim_scores(sim_scores):
    # Set raw score, which is used by Orphan placement
    sim_scores['raw_score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
    # Distribute final scores_components between 0-1.
    sim_scores['psi'] = (sim_scores['psi'] - sim_scores['psi'].min()) / \
                        (sim_scores['psi'].max() - sim_scores['psi'].min())

    sim_scores['subsmat'] = (sim_scores['subsmat'] - sim_scores['subsmat'].min()) / \
                            (sim_scores['subsmat'].max() - sim_scores['subsmat'].min())

    # ToDo: Experiment testing these magic number weights...
    sim_scores['score'] = (sim_scores['psi'] * 0.3) + (sim_scores['subsmat'] * 0.7)
    return sim_scores


def mc_create_sequence_hmms(record, args):
    outdir = args[0]
    # First check to see if they already exist
    hmm_dir = join(outdir, "hmm")
    if not os.path.isfile(join(hmm_dir, "%s.hmm" % record.id)):
        align = Alb.AlignBuddy(record.format("fasta"))
        align = Alb.generate_hmm(align, HMMBUILD)
        with open(join(hmm_dir, "%s.hmm" % record.id), "w") as _ofile:
            _ofile.write(align.alignments[0].hmm)
    return


class WorkerJob(object):
    def __init__(self, seqbuddy, sql_broker):
        self.seqbuddy = seqbuddy
        self.seq_ids = sorted([rec.id for rec in self.seqbuddy.records])
        self.seq_id_hash = hlp.md5_hash(", ".join(self.seq_ids))
        # Job id can't be the same as seq_id_hash, because different alignment options may be involved
        job_id = "".join([str(x) for x in [self.seq_id_hash, GAP_OPEN, GAP_EXTEND, ALIGNMETHOD, ALIGNPARAMS, TRIMAL]])
        self.job_id = hlp.md5_hash(job_id)
        self.sql_broker = sql_broker
        self.heartbeat = HeartBeat(HEARTBEAT_DB, MASTER_PULSE)
        self.running = False

    def run(self):
        self.heartbeat.start()
        with hlp.ExclusiveConnect(HEARTBEAT_DB) as cursor:
            workers = cursor.execute("SELECT * FROM heartbeat "
                                     "WHERE thread_type='worker' "
                                     "AND pulse>%s" % (time.time() - MAX_WORKER_WAIT - cursor.lag)).fetchone()
        result = False
        if workers:
            self.running = True
            self.queue_job()
            while self.running:
                db_result = self.pull_from_db()
                if db_result:
                    result = db_result
                    break

                if self.check_finished():
                    result = self.process_finished()
                    if result:
                        result = result
                        break

                elif random() > 0.75:  # Occasionally check to make sure the job is still active
                    if not self.check_if_active():
                        break

                with hlp.ExclusiveConnect(HEARTBEAT_DB) as cursor:
                    active_threads = cursor.execute("SELECT COUNT(*) FROM heartbeat "
                                                    "WHERE thread_type='master'").fetchone()[0]
                pause_time = active_threads * 0.5
                time.sleep(pause_time)
        self.heartbeat.end()
        return result

    def queue_job(self):
        with hlp.ExclusiveConnect(WORKER_DB) as cursor:
            # Make sure the job hasn't already been submitted or completed)
            complete_check = cursor.execute("SELECT hash FROM complete WHERE hash=?", (self.job_id,)).fetchone()
            queue_check = cursor.execute("SELECT hash FROM queue WHERE hash=?", (self.job_id,)).fetchone()
            processing_check = cursor.execute("SELECT hash FROM processing WHERE hash=?", (self.job_id,)).fetchone()

            if not complete_check and not queue_check and not processing_check:
                if not os.path.isfile("%s/%s.seqs" % (WORKER_OUT, self.job_id)):
                    self.seqbuddy.write("%s/%s.seqs" % (WORKER_OUT, self.job_id), out_format="fasta")
                cursor.execute("INSERT INTO queue "
                               "(hash, psi_pred_dir, align_m, align_p, trimal, gap_open, gap_extend) "
                               "VALUES (?, ?, ?, ?, ?, ?, ?)",
                               (self.job_id, PSIPREDDIR, ALIGNMETHOD, ALIGNPARAMS, " ".join([str(x) for x in TRIMAL]),
                                GAP_OPEN, GAP_EXTEND,))

            if not cursor.execute("SELECT * FROM waiting WHERE hash=? AND master_id=?",
                                  (self.job_id, self.heartbeat.id,)).fetchone():
                cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", (self.job_id, self.heartbeat.id,))
        return

    def pull_from_db(self):
        # Remember that job_id and seq_id_hash are different! job_id includes details about alignment
        query = self.sql_broker.query("SELECT graph, alignment FROM data_table WHERE hash=?", (self.seq_id_hash,))
        if not query or len(query[0]) != 2:
            return False
        sim_scores, alignment = query[0]
        sim_scores = pd.read_csv(StringIO(sim_scores), index_col=False, header=None)
        sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
        with hlp.ExclusiveConnect(WORKER_DB) as cursor:
            cursor.execute("DELETE FROM waiting WHERE hash=? AND master_id=?", (self.job_id,
                                                                                self.heartbeat.id,))
        self.heartbeat.end()
        return sim_scores, Alb.AlignBuddy(alignment, in_format="fasta")

    def check_finished(self):
        with hlp.ExclusiveConnect(WORKER_DB) as cursor:
            if not cursor.execute("SELECT * FROM complete WHERE hash=?", (self.job_id,)).fetchone():
                return False

            cursor.execute("DELETE FROM complete WHERE hash=?", (self.job_id,))
            cursor.execute("INSERT INTO proc_comp (hash, master_id) VALUES (?, ?)",
                           (self.job_id, self.heartbeat.id))
            cursor.execute("DELETE FROM waiting WHERE hash=? AND master_id=?", (self.job_id, self.heartbeat.id,))
        return True

    def process_finished(self):
        location = join(WORKER_OUT, self.job_id)
        alignment = Alb.AlignBuddy("%s.aln" % location, in_format="fasta")
        sim_scores = pd.read_csv("%s.graph" % location, index_col=False, header=None)
        sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
        cluster2database(Cluster(self.seq_ids, sim_scores), self.sql_broker, alignment)

        for del_file in [".aln", ".graph", ".seqs"]:
            os.remove(join(WORKER_OUT, "%s%s" % (self.job_id, del_file)))
        with hlp.ExclusiveConnect(WORKER_DB) as cursor:
            cursor.execute("DELETE FROM proc_comp WHERE hash=?", (self.job_id,))
        return sim_scores, alignment

    def check_if_active(self):
        with hlp.ExclusiveConnect(HEARTBEAT_DB, priority=True) as hb_cursor:
            min_pulse = time.time() - MAX_WORKER_WAIT - hb_cursor.lag
            worker_heartbeat_check = hb_cursor.execute("SELECT * FROM heartbeat "
                                                       "WHERE thread_type='worker' AND pulse>%s"
                                                       % min_pulse).fetchall()
            master_heartbeat_check = hb_cursor.execute("SELECT * FROM heartbeat "
                                                       "WHERE thread_type='master' AND pulse>%s"
                                                       % min_pulse).fetchall()

        with hlp.ExclusiveConnect(WORKER_DB) as cursor:
            if not worker_heartbeat_check:
                # No workers are still around. Time to clean up and move on serially
                cursor.execute("DELETE FROM queue WHERE hash=?", (self.job_id,))
                cursor.execute("DELETE FROM waiting WHERE master_id=?", (self.heartbeat.id,))
                try:
                    os.remove(join(WORKER_OUT, "%s.seqs" % self.job_id))
                except FileNotFoundError:
                    pass
                return False

            queue_check = cursor.execute("SELECT * FROM queue WHERE hash=?", (self.job_id,)).fetchone()
            processing_check = cursor.execute("SELECT * FROM processing WHERE hash=?",
                                              (self.job_id,)).fetchone()
            complete_check = cursor.execute("SELECT * FROM complete WHERE hash=?", (self.job_id,)).fetchone()
            comp_proc_check = cursor.execute("SELECT * FROM proc_comp WHERE hash=?", (self.job_id,)).fetchone()

        if processing_check:
            # Make sure the respective worker is still alive to finish processing this job
            workers = [(thread_id, pulse) for thread_id, thread_type, pulse in worker_heartbeat_check]
            workers = OrderedDict(workers)
            worker_id = processing_check[1]
            if worker_id not in workers:
                # Doesn't look like it... Might as well queue it back up
                self.queue_job()

        elif comp_proc_check:
            # Make sure respective master is still alive to finish processing this job
            masters = [(thread_id, pulse) for thread_id, thread_type, pulse in master_heartbeat_check]
            masters = OrderedDict(masters)
            master_id = comp_proc_check[1]
            if master_id not in masters:
                # Nope... Put it back into complete and pick it up on the next loop
                with hlp.ExclusiveConnect(WORKER_DB) as cursor:
                    cursor.execute("INSERT INTO complete (hash) VALUES (?)", (self.job_id,))
                    cursor.execute("DELETE FROM proc_comp WHERE hash=?", (self.job_id,))
        elif queue_check:
            pass

        elif complete_check:
            pass

        else:
            # Check again if job is in main DB, it is only deleted from work_db after it has been added to main DB
            if not self.pull_from_db():
                # Ummmm... Where'd the job go??? Queue it back up because it really seems to have vanished
                print("%s vanished!" % self.job_id)
                self.queue_job()
        return True

# ################ END SCORING FUNCTIONS ################ #


# #########  MCL stuff  ########## #
def mcmcmc_mcl(args, params):
    """
    Function passed to mcmcmcm.MCMCMC that will execute MCL and return the final cluster scores
    :param args: Sample values to run MCL with and a random seed [inflation, gq, r_seed]
    :param params: List of parameters (see below for unpacking assignment)
    :return:
    """
    inflation, gq, r_seed = args
    exter_tmp_dir, seqbuddy, parent_cluster, taxa_sep, \
        sql_broker, psi_pred_ss2, progress, expect_num_results = params
    rand_gen = Random(r_seed)
    mcl_obj = hlp.MarkovClustering(parent_cluster.sim_scores, inflation=inflation, edge_sim_threshold=gq)
    mcl_obj.run()
    progress.update('mcl_runs', 1)
    clusters = mcl_obj.clusters
    # Order the clusters so the big jobs are queued up front.
    clusters = sorted(clusters, key=lambda x: len(x), reverse=True)
    score = 0

    child_list = OrderedDict()
    for indx, cluster_ids in enumerate(clusters):
        sb_copy = Sb.make_copy(seqbuddy)
        sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in cluster_ids]))
        # Queue jobs if appropriate
        if WORKER_DB and os.path.isfile(WORKER_DB) and len(cluster_ids) >= MIN_SIZE_TO_WORKER:
            p = Process(target=mc_create_all_by_all_scores, args=(sb_copy, [psi_pred_ss2, sql_broker]))
            # Need to slow the start down a little, otherwise it can run into concurrency issues
            time.sleep(0.1)
            p.start()
            seq_ids = sorted([rec.id for rec in sb_copy.records])
            seq_id_hash = hlp.md5_hash(", ".join(seq_ids))
            child_list[seq_id_hash] = [p, indx, cluster_ids]
        else:
            sim_scores, alb_obj = retrieve_all_by_all_scores(sb_copy, psi_pred_ss2, sql_broker, quiet=True)
            cluster = Cluster(cluster_ids, sim_scores, parent=parent_cluster, taxa_sep=taxa_sep,
                              r_seed=rand_gen.randint(1, 999999999999999))
            clusters[indx] = cluster
            score += cluster.score()

    # wait for remaining processes to complete
    while len(child_list) > 0:
        for _name, args in child_list.items():
            child, indx, cluster_ids = args
            if child.is_alive():
                continue
            else:
                query = sql_broker.query("SELECT graph, alignment FROM data_table WHERE hash=?", (_name,))
                if query and len(query[0]) == 2:
                    sim_scores, alignment = query[0]
                    alignment = Alb.AlignBuddy(alignment, in_format="fasta")
                    if len(alignment.records()) == 1:
                        sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
                    else:
                        sim_scores = pd.read_csv(StringIO(sim_scores), index_col=False, header=None)
                        sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]

                    cluster = Cluster(cluster_ids, sim_scores, parent=parent_cluster, taxa_sep=taxa_sep,
                                      r_seed=rand_gen.randint(1, 999999999999999))
                    clusters[indx] = cluster
                    score += cluster.score()
                    del child_list[_name]
                    break

    with LOCK:
        with open(join(exter_tmp_dir, "max.txt"), "r") as ifile:
            results = ifile.readlines()
            results = [result.strip() for result in results]
            results.append(",".join([cluster.seq_id_hash for cluster in clusters]))
            results = sorted(results)

        with open(join(exter_tmp_dir, "max.txt"), "w") as ofile:
            ofile.write("\n".join(results))

    if len(results) == expect_num_results:
        best_score = None
        best_clusters = []  # Hopefully just find a single best set of cluster, but could be more
        for clusters in results:
            score_sum = 0
            cluster_ids = []
            for cluster in clusters.split(","):
                sql_query = sql_broker.query("SELECT seq_ids, graph FROM data_table WHERE hash=?", (cluster,))
                seq_ids = sql_query[0][0].split(", ")
                cluster_ids.append(sql_query[0][0])
                if len(seq_ids) == 1:
                    sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
                else:
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
        with LOCK:
            with open(join(exter_tmp_dir, "best_group"), "w") as ofile:
                ofile.write('\n'.join(best_clusters))
            open(join(exter_tmp_dir, "max.txt"), "w").close()
    elif len(results) > expect_num_results:  # This should never be able to happen
        raise ValueError("More results written to max.txt than expect_num_results")
    return score


def parse_mcl_clusters(path):
    with open(path, "r") as ifile:
        clusters = ifile.readlines()
    clusters = [cluster.strip().split("\t") for cluster in clusters]
    return clusters


def write_mcl_clusters(clusters, path):
    clusters_strings = ["\t".join(sorted(cluster.seq_ids)) for cluster in clusters]
    with open(path, "w") as ofile:
        ofile.write("\n".join(clusters_strings))
    return


'''
# #########  Old orphan placement  (deprecated) ########## #
class Orphans(object):  # Deprecated
    def __init__(self, seqbuddy, clusters, sql_broker, psi_pred_ss2, outdir, min_clust_size=4, quiet=False):
        """
        Organizes all of the orphan prediction/folding logic into a class
        :param seqbuddy: This should include all of the sequences in the entire population being tests
        :param clusters: List of all cluster objects the sequences have currently been grouped into
        :param sql_broker: Open SQLiteBroker object
        :param psi_pred_ss2: OrderedDict of all PSI Pred graph paths (seq_id: df)
        :param outdir: Primary location of RD-MCL run
        :param min_clust_size: Specify the smallest group that can be classified as an actual 'cluster'
        :param quiet: Suppress any terminal output
        """
        self.seqbuddy = seqbuddy
        self.clusters = clusters
        self.sql_broker = sql_broker
        self.psi_pred_ss2 = psi_pred_ss2
        self.outdir = outdir
        self.min_clust_size = min_clust_size
        self.quiet = quiet
        self.num_orphans = 0
        self.small_clusters = OrderedDict()
        self.large_clusters = OrderedDict()
        self._separate_large_small()
        self.printer = br.DynamicPrint(quiet=self.quiet)
        self.tmp_file = br.TempFile()
        self.tmp_dir = br.TempDir()
        self.tmp_dir.subfile("seqs.fa")
        self.seqbuddy.write(self.tmp_dir.subfiles[0], "fasta")

        # Create HMM forward-score dataframe from initial input --> p(seq|hmm)
        self.hmm_fwd_scores = pd.DataFrame(columns=["group", "rec_id", "fwd_raw"])
        for group_name in self.large_clusters:
            hmm_path = join(self.outdir, "hmm", group_name)
            fwdback_output = Popen("generic_fwdback_example %s %s" % (hmm_path, self.tmp_dir.subfiles[0]),
                                   shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].decode()
            fwd_scores_df = pd.read_csv(StringIO(fwdback_output), delim_whitespace=True,
                                        header=None, comment="#", index_col=False)
            fwd_scores_df.columns = ["rec_id", "fwd_raw", "back_raw", "fwd_bits", "back_bits"]
            fwd_scores_df["group"] = group_name
            self.hmm_fwd_scores = self.hmm_fwd_scores.append(fwd_scores_df.loc[:, ["group", "rec_id", "fwd_raw"]],
                                                             ignore_index=True)

        # Calculate all-by-all matrix of correlation coefficients on fwd scores among all sequences
        self.rsquare_vals_df = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
        sub_recs = list(seqbuddy.records[1:])
        for rec1 in seqbuddy.records:
            for rec2 in sub_recs:
                fwd1 = self.hmm_fwd_scores.loc[self.hmm_fwd_scores.rec_id == rec1.id].sort_values(by="group").fwd_raw
                fwd2 = self.hmm_fwd_scores.loc[self.hmm_fwd_scores.rec_id == rec2.id].sort_values(by="group").fwd_raw
                corr = scipy.stats.pearsonr(fwd1, fwd2)
                comparison = pd.DataFrame(data=[[rec1.id, rec2.id, corr[0]**2]],
                                          columns=["rec_id1", "rec_id2", "r_square"])
                self.rsquare_vals_df = self.rsquare_vals_df.append(comparison, ignore_index=True)
            sub_recs = sub_recs[1:]

        self.lrg_cluster_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square", "group"])
        for group_name, clust in self.large_clusters.items():
            df = self.rsquare_vals_df.loc[(self.rsquare_vals_df["rec_id1"].isin(clust.seq_ids)) &
                                          (self.rsquare_vals_df["rec_id2"].isin(clust.seq_ids))].copy()
            df["group"] = group_name
            self.lrg_cluster_rsquares = self.lrg_cluster_rsquares.append(df)

    def _separate_large_small(self):
        for cluster in self.clusters:
            if len(cluster.seq_ids) < self.min_clust_size:
                self.small_clusters[cluster.name()] = cluster
            else:
                self.large_clusters[cluster.name()] = cluster

    def mc_temp_merge_clusters(self, args):
        clust1, clust2 = args
        self.temp_merge_clusters(clust1, clust2)
        return

    def temp_merge_clusters(self, clust1, clust2):
        # Read or create an alignment containing the sequences in the clusters being merged
        seq_ids = sorted(clust1.seq_ids + clust2.seq_ids)
        seqbuddy = Sb.make_copy(self.seqbuddy)
        regex = "^%s$" % "$|^".join(seq_ids)
        Sb.pull_recs(seqbuddy, regex)

        sim_scores, alb_obj = retrieve_all_by_all_scores(seqbuddy, self.psi_pred_ss2, self.sql_broker, quiet=True)
        cluster2database(Cluster(seq_ids, sim_scores), self.sql_broker, alb_obj)
        return sim_scores

    def _check_orphan(self, small_cluster):
        # First prepare the data for each large cluster so they can be compared
        log_output = "%s\n" % small_cluster.seq_ids
        log_output += "Min sim score needed: %s\n" % self.lrg_cluster_rsquares.r_square.min()
        data_dict = OrderedDict()
        for group_name, large_cluster in self.large_clusters.items():
            log_output += "\t%s\n" % group_name
            log_output += "\t%s\n" % large_cluster.seq_ids

            lrg2lrg_data = self.rsquare_vals_df.loc[(self.rsquare_vals_df["rec_id1"].isin(large_cluster.seq_ids))
                                                    & (self.rsquare_vals_df["rec_id2"].isin(large_cluster.seq_ids))]

            log_output += "\tAmong large:\t%s - %s\n" % (round(lrg2lrg_data.r_square.min(), 4),
                                                         round(lrg2lrg_data.r_square.max(), 4))

            # Pull out the similarity scores between just the large and small cluster sequences
            lrg2sml_data = self.rsquare_vals_df.loc[(self.rsquare_vals_df["rec_id1"].isin(small_cluster.seq_ids))
                                                    | (self.rsquare_vals_df["rec_id2"].isin(small_cluster.seq_ids))]

            lrg2sml_data = lrg2sml_data.loc[(lrg2sml_data['rec_id1'].isin(large_cluster.seq_ids))
                                            | (lrg2sml_data['rec_id2'].isin(large_cluster.seq_ids))]

            log_output += "\tLarge 2 small:\t%s - %s\n\n" % (round(lrg2sml_data.r_square.min(), 4),
                                                             round(lrg2sml_data.r_square.max(), 4))

            # Confirm that the orphans are sufficiently similar to large group to warrant consideration by ensuring
            # that at least one similarity score is greater than the lowest score with all large groups.
            # Also, convert group_data to a numpy array so sm.stats can read it
            if lrg2sml_data.r_square.max() >= self.lrg_cluster_rsquares.r_square.min():
                data_dict[group_name] = (True, np.array(lrg2sml_data.r_square))
            else:
                data_dict[group_name] = (False, np.array(lrg2sml_data.r_square))
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
            mean_diff = abs(np.mean(data_dict[max_ave_name][1]) - np.mean(self.lrg_cluster_rsquares.r_square))
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
            mean_diff = abs(np.mean(data_dict[max_ave_name][1]) - np.mean(self.lrg_cluster_rsquares.r_square))
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

    def place_orphans(self):
        if not self.small_clusters or not self.large_clusters or len(self.large_clusters) == 1:
            return
        best_cluster = OrderedDict([("small_name", None), ("large_name", None), ("meandiff", 0)])
        fostered_orphans = 0
        starting_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])

        # Loop through repeatedly, adding a single small cluster at a time and updating the large clusters on each loop
        while True and self.large_clusters:
            remaining_orphans = sum([len(sm_clust.seq_ids) for group, sm_clust in self.small_clusters.items()])
            self.printer.write("%s of %s orphans remain" % (remaining_orphans, starting_orphans))
            small_clusters = [clust for clust_name, clust in self.small_clusters.items()]

            for clust in small_clusters:
                foster_score = self._check_orphan(clust)
                if foster_score:
                    large_name, mean_diff = foster_score
                    if mean_diff > best_cluster["meandiff"]:
                        best_cluster = OrderedDict([("small_name", clust.name()), ("large_name", large_name),
                                                    ("meandiff", mean_diff)])

            if best_cluster["small_name"]:
                small_name = best_cluster["small_name"]
                large_name = best_cluster["large_name"]
                large_cluster = self.large_clusters[large_name]
                small_cluster = self.small_clusters[small_name]
                large_cluster.reset_seq_ids(small_cluster.seq_ids + large_cluster.seq_ids)

                seqs = Sb.pull_recs(Sb.make_copy(self.seqbuddy), "^%s$" % "$|^".join(large_cluster.seq_ids))
                sim_scores, alb_obj = retrieve_all_by_all_scores(seqs, self.psi_pred_ss2,
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
'''


class Seqs2Clusters(object):
    def __init__(self, clusters, min_clust_size, seqbuddy, outdir):
        self.clusters = clusters
        self.min_clust_size = min_clust_size
        self.seqbuddy = seqbuddy
        self.outdir = outdir
        self.tmp_file = br.TempFile()
        self.tmp_dir = br.TempDir()
        self.small_clusters = OrderedDict()
        self.large_clusters = OrderedDict()
        self._separate_large_small()

    def _separate_large_small(self):
        for cluster in self.clusters:
            if len(cluster.seq_ids) < self.min_clust_size:
                self.small_clusters[cluster.name()] = cluster
            else:
                self.large_clusters[cluster.name()] = cluster
        return

    def _mc_build_cluster_nulls(self, clust, args):
        rsquare_vals_df, global_null_file, cluster_nulls_file, out_of_cluster_file, temp_log_output = args

        log_output = "%s\n" % clust.name()
        clust_null_df = rsquare_vals_df.loc[(rsquare_vals_df["seq1"].isin(clust.seq_ids)) &
                                            (rsquare_vals_df["seq2"].isin(clust.seq_ids))].copy()
        global_null_output = clust_null_df.to_csv(path_or_buf=None, header=None, index=False, index_label=False)
        cluster_nulls_output = ""
        clust_null_df = clust_null_df.loc[clust_null_df.seq1 != clust_null_df.seq2].reset_index(drop=True)

        out_of_cluster_output = ""
        if len(clust_null_df) > 1:
            cluster_nulls_output += '"%s":{"mu":%s,"sigma":%s},' % (clust.name(), hlp.mean(clust_null_df.r_square),
                                                                    hlp.std(clust_null_df.r_square))
            log_output += "\tN: %s\n" % len(clust_null_df)
            log_output += "\tMean: %s\n" % hlp.mean(clust_null_df.r_square)
            log_output += "\tStd: %s\n\n" % hlp.std(clust_null_df.r_square)
        else:
            log_output += "\tN: 1\n"
            log_output += "\tMean: Null\n"
            log_output += "\tStd: Null\n\n"

        if clust.name() in self.large_clusters:
            clust_null_df = rsquare_vals_df.loc[((rsquare_vals_df["seq1"].isin(clust.seq_ids)) &
                                                 (-rsquare_vals_df["seq2"].isin(clust.seq_ids))) |
                                                ((-rsquare_vals_df["seq1"].isin(clust.seq_ids)) &
                                                 (rsquare_vals_df["seq2"].isin(clust.seq_ids)))].copy()
            out_of_cluster_output += clust_null_df.to_csv(path_or_buf=None, header=None, index=False, index_label=False)

        with LOCK:
            with open(global_null_file, "a") as ofile:
                ofile.write(global_null_output)
            with open(cluster_nulls_file, "a") as ofile:
                ofile.write(cluster_nulls_output)
            with open(out_of_cluster_file, "a") as ofile:
                ofile.write(out_of_cluster_output)
            with open(temp_log_output, "a") as ofile:
                ofile.write(log_output)
        return

    def _mc_build_seq2group(self, seq_id, args):
        rsquare_vals_df, seq2group_dists_file, orig_clusters_file, temp_log_output = args
        seq2group_dists = '"%s":{' % seq_id
        orig_clusters = '"%s":' % seq_id
        log_output = "%s\n" % seq_id
        best_hit = None
        for clust in self.clusters:
            if seq_id in clust.seq_ids:
                orig_clusters += '"%s",' % clust.name()
            log_output += "\t" + clust.name() + ": "
            seq_ids = set([i for i in clust.seq_ids if i != seq_id])

            # Pull out dataframe rows where current sequence is paired with a sequence from the current cluster
            df = rsquare_vals_df.loc[((rsquare_vals_df["seq1"].isin(seq_ids)) |
                                      (rsquare_vals_df["seq2"].isin(seq_ids))) &
                                     ((rsquare_vals_df["seq1"] == seq_id) |
                                      (rsquare_vals_df["seq2"] == seq_id)) &
                                     (rsquare_vals_df["seq1"] != rsquare_vals_df["seq2"])].copy()
            if df.empty:
                test_mean = 0
            else:
                test_mean = hlp.mean(df.r_square)

            seq2group_dists += '"%s":%s,' % (clust.name(), test_mean)
            log_output += "%s\n" % test_mean
            if not best_hit:
                best_hit = (clust.name(), test_mean)
            else:
                if test_mean > best_hit[1]:
                    best_hit = (clust.name(), test_mean)

        log_output += "\tBest: %s\n" % best_hit[0]
        seq2group_dists = seq2group_dists.strip(",") + "},"
        with LOCK:
            with open(seq2group_dists_file, "a") as ofile:
                ofile.write("%s" % seq2group_dists)
            with open(orig_clusters_file, "a") as ofile:
                ofile.write(orig_clusters)
            with open(temp_log_output, "a") as ofile:
                ofile.write(log_output)
        return

    @staticmethod
    def _create_truncnorm(mu, sigma, lower=0, upper=1):
        sigma = sigma if sigma > 0.001 else 0.001  # This prevents unrealistically small differences and DivBy0 errors
        dist = scipy.stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
        return dist

    @staticmethod
    def get_sort_order(seq2group_dists, placed, new_group, orig_clusts):
        """
        Arrange seq_ids by cluster size (small to big), then by R² within each size bin (high to low)
        :param seq2group_dists: {seq_id: {clust_name: R² mean, ...}, ...}
        :param placed: [seq_id, ...]
        :param new_group: {clust_name: {seq_id, ...}, ...}
        :param orig_clusts: {seq_id: clust_name, ...}
        :return: [[(seq_id, R² mean), ...], ...]
        """
        sort_vals = [(seq_id, max([val for group_id, val in group_vals.items()]), len(new_group[orig_clusts[seq_id]]))
                     for seq_id, group_vals in seq2group_dists.items()]
        sort_vals = [(seq_id, val, len_group) for seq_id, val, len_group in sorted(sort_vals, key=lambda i: i[2])
                     if seq_id not in placed]
        sort_order = [[]]
        current_group_len = 0
        valve = br.SafetyValve()
        for seq_id, val, len_group in sort_vals:
            while current_group_len < len_group:
                valve.step()
                if sort_order[-1]:
                    sort_order.append([])
                current_group_len += 1
            else:
                sort_order[-1].append((seq_id, val))
            sort_order[-1] = sorted(sort_order[-1], key=lambda i: i[1], reverse=True)
        return sort_order

    def place_seqs_in_clusts(self):
        """
        Create null distributions for each existing group, then check all sequences against all groups
        :return:
        """
        # Calculate HMMs from every sequence, calculate the forward probability for every sequences against every HMM,
        # and then calculate correlation coefficients between every pair of sequences (i.e., make an all-by-all matrix)
        log_output = br.TempFile()
        log_output.write("# HMM forward scores all-by-all R² dataframe #\n")
        fwd_scores = FwdScoreCorrelations(self.seqbuddy, self.outdir)
        rsquare_vals_df = fwd_scores.rsquare_vals_df
        '''
        rsquares_df_path = join(self.outdir, "hmm", "rsquares_matrix.csv")
        if os.path.isfile(rsquares_df_path):
            log_output.write("\tRead from %s\n\n" % rsquares_df_path)
            rsquare_vals_df = pd.read_csv(rsquares_df_path)
        else:
            log_output.write("\tCalculated new and written to %s\n\n" % rsquares_df_path)
            rsquare_vals_df = self.create_fwd_score_rsquared_matrix()
            rsquare_vals_df.to_csv(rsquares_df_path, index=False)

        hmm_fwd_scores_path = join(self.outdir, "hmm", "hmm_fwd_scores.csv")
        if not os.path.isfile(hmm_fwd_scores_path):
            hmm_fwd_scores = self.create_hmm_fwd_score_df()
            hmm_fwd_scores.to_csv(hmm_fwd_scores_path, index=False)
        '''
        valve = br.SafetyValve(round(len(self.clusters[0].get_base_cluster().seq_ids) * 1.2))
        placed = []

        while True:
            try:
                valve.step()
            except RuntimeError:
                log_output.write("\n\nERROR: place_seqs_in_clusts blew its safety valve; abort further placement.\n")
                break

            # Create a copy of the input clusters for checking later
            copy_clusters = deepcopy(self.clusters)

            if len(placed) == len(self.clusters[0].get_base_cluster()):
                placed = []  # Reset the `placed` list and loop through again, making sure nothing has changed
                breakout = True
            else:
                breakout = False

            log_output.write("# Current clusters #\n")
            for clust in self.clusters:
                log_output.write("%s\n%s\n\n" % (clust.name(), clust.seq_ids))

            # Create a distribution from all R² values between records in each pre-called large cluster
            global_null_file = br.TempFile()
            global_null_file.write("seq1,seq2,r_square\n")

            out_of_cluster_file = br.TempFile()
            out_of_cluster_file.write("seq1,seq2,r_square\n")

            # Also calculate sub-distributions from within each pre-called cluster, creating a cluster null
            log_output.write("# Cluster R² distribution statistics #\n")
            cluster_nulls_file = br.TempFile()
            cluster_nulls_file.write("{")

            temp_log_output = br.TempFile()

            args = [rsquare_vals_df, global_null_file.path, cluster_nulls_file.path,
                    out_of_cluster_file.path, temp_log_output.path]
            br.run_multicore_function(self.clusters, self._mc_build_cluster_nulls, args, max_processes=CPUS, quiet=True)

            global_null_df = pd.read_csv(global_null_file.path)
            try:
                cluster_nulls = json.loads(cluster_nulls_file.read().strip(",") + "}")
            except json.JSONDecodeError as err:
                print(cluster_nulls_file.read())
                raise err
            cluster_nulls = {clust_name: self._create_truncnorm(dist["mu"], dist["sigma"])
                             for clust_name, dist in cluster_nulls.items()}
            out_of_cluster_df = pd.read_csv(out_of_cluster_file.path)
            log_output.write(temp_log_output.read())
            temp_log_output.clear()

            # H₀: R² between two sequences is drawn from the same distribution as pre-called cluster pairs
            global_null_df = global_null_df.loc[global_null_df.seq1 != global_null_df.seq2].reset_index(drop=True)
            log_output.write("# Global null R² distribution statistics #\n")
            log_output.write("\tN: %s\n" % len(global_null_df))
            log_output.write("\tMean: %s\n" % hlp.mean(global_null_df.r_square))
            log_output.write("\tStd: %s\n\n" % hlp.std(global_null_df.r_square))

            out_of_cluster_df = out_of_cluster_df.reset_index(drop=True)
            log_output.write("# Out of cluster R² distribution statistics #\n")
            log_output.write("\tN: %s\n" % len(out_of_cluster_df))
            log_output.write("\tMean: %s\n" % hlp.mean(out_of_cluster_df.r_square))
            log_output.write("\tStd: %s\n\n" % hlp.std(out_of_cluster_df.r_square))

            # Create Gaussian distribution object that is clipped off below 0 and above 1
            null_dist = self._create_truncnorm(hlp.mean(global_null_df.r_square),
                                               hlp.std(global_null_df.r_square))
            out_of_cluster_dist = self._create_truncnorm(hlp.mean(out_of_cluster_df.r_square),
                                                         hlp.std(out_of_cluster_df.r_square))

            # Calculate how well every sequence fits with every group
            log_output.write("# Group placements #\n")
            seq2group_dists_file = br.TempFile()
            seq2group_dists_file.write("{")

            orig_clusters_file = br.TempFile()
            orig_clusters_file.write("{")

            args = [rsquare_vals_df, seq2group_dists_file.path, orig_clusters_file.path, temp_log_output.path]
            seq_ids = self.clusters[0].get_base_cluster().seq_ids
            br.run_multicore_function(seq_ids, self._mc_build_seq2group, args, quiet=True, max_processes=CPUS)

            seq2group_dists = json.loads(seq2group_dists_file.read().strip(",") + "}")
            orig_clusters = json.loads(orig_clusters_file.read().strip(",") + "}")
            log_output.write(temp_log_output.read())

            new_groups = OrderedDict([(clust.name(), copy(clust.seq_ids)) for clust in self.clusters])
            new_groups["orphans"] = set()

            # Set a static order that sequences will be placed. Smallest clusters and largest average r_square first.
            sort_order = self.get_sort_order(seq2group_dists, placed, new_groups, orig_clusters)
            log_output.write("\n\n# Sort order #\n%s\n\n" % sort_order)
            sort_order = [seq_id for group in sort_order for seq_id, val in group]

            log_output.write("# Placement #\n")
            for seq_id in sort_order:
                highest = [(group_id, val) for group_id, val in seq2group_dists[seq_id].items()]
                group_id, val = sorted(highest, key=lambda i: i[1], reverse=True)[0]
                log_output.write("\t%s\t%s\t%s\n" % (seq_id, group_id, val))
                log_output.write("\tOut of Clust = %s\n" % out_of_cluster_dist.pdf(val))
                log_output.write("\tGlobal null = %s\n " % null_dist.pdf(val))
                log_output.write("\tCluster null = NONE\n" if group_id not in cluster_nulls
                                 else "\tCluster null = %s\n" % cluster_nulls[group_id].pdf(val))

                # If the 'highest' cluster is not the original cluster and is at least 90% smaller than the original
                # cluster, then check the original cluster first
                if group_id != orig_clusters[seq_id] \
                        and orig_clusters[seq_id] in cluster_nulls \
                        and len(new_groups[group_id]) <= len(new_groups[orig_clusters[seq_id]]) * 0.1:
                    orig_val = seq2group_dists[seq_id][orig_clusters[seq_id]]
                    if cluster_nulls[orig_clusters[seq_id]].pdf(orig_val) >= out_of_cluster_dist.pdf(orig_val) \
                            or null_dist.pdf(orig_val) >= out_of_cluster_dist.pdf(orig_val):
                        group_id, val = orig_clusters[seq_id], orig_val

                # Place sequence in new_groups and update seq2group_dists
                orphaned = False
                if group_id not in cluster_nulls and null_dist.pdf(val) < out_of_cluster_dist.pdf(val):
                    orphaned = True
                elif group_id in cluster_nulls \
                        and cluster_nulls[group_id].pdf(val) < out_of_cluster_dist.pdf(val)\
                        and null_dist.pdf(val) < out_of_cluster_dist.pdf(val):
                    orphaned = True

                if orphaned:
                    log_output.write("\tOrphaned\n")
                    new_groups["orphans"].add(seq_id)

                placed.append(seq_id)
                if seq_id in new_groups["orphans"] or seq_id not in new_groups[group_id]:  # The Sequence has moved
                    # Add
                    if seq_id not in new_groups["orphans"]:
                        log_output.write("\tMoving from %s to %s\n\n" % (orig_clusters[seq_id], group_id))
                        new_groups[group_id].add(seq_id)

                    # Remove
                    if len(new_groups[orig_clusters[seq_id]]) == 1 and seq_id in new_groups["orphans"]:
                        # If the sequence started orphaned, don't change its group ID
                        log_output.write("\tRetaining original orphan %s\n\n" % orig_clusters[seq_id])
                        new_groups["orphans"].remove(seq_id)

                    else:
                        orig_clust = orig_clusters[seq_id]
                        log_output.write("\tCreating new group\n\n" if seq_id in new_groups["orphans"] else "")
                        new_groups[orig_clust].remove(seq_id)

                        if not new_groups[orig_clust]:  # When removing a sequence leaves the cluster empty
                            del new_groups[orig_clust]
                            for subseq_id, group_vals in copy(seq2group_dists).items():
                                if orig_clust in group_vals:
                                    del seq2group_dists[subseq_id][orig_clust]

                        else:
                            # Update against old cluster
                            seq_ids = new_groups[orig_clust]
                            for si in sort_order:
                                df = rsquare_vals_df.loc[((rsquare_vals_df["seq1"].isin(seq_ids)) |
                                                          (rsquare_vals_df["seq2"].isin(seq_ids))) &
                                                         ((rsquare_vals_df["seq1"] == si) |
                                                          (rsquare_vals_df["seq2"] == si)) &
                                                         (rsquare_vals_df["seq1"] != rsquare_vals_df["seq2"])]
                                df = df.copy()

                                if df.empty:
                                    test_mean = 0
                                else:
                                    test_mean = hlp.mean(df.r_square)
                                seq2group_dists[si][orig_clust] = test_mean

                        # Update against new cluster
                        seq_ids = new_groups[group_id]
                        for si in sort_order:
                            df = rsquare_vals_df.loc[((rsquare_vals_df["seq1"].isin(seq_ids)) |
                                                      (rsquare_vals_df["seq2"].isin(seq_ids))) &
                                                     ((rsquare_vals_df["seq1"] == si) |
                                                      (rsquare_vals_df["seq2"] == si)) &
                                                     (rsquare_vals_df["seq1"] != rsquare_vals_df["seq2"])].copy()
                            if df.empty:
                                test_mean = 0
                            else:
                                test_mean = hlp.mean(df.r_square)
                            seq2group_dists[si][group_id] = test_mean
                        breakout = False
                        break
                else:  # The sequence was returned to its original group
                    log_output.write("\tReturned to group %s\n\n" % group_id)

            for clust in self.clusters:
                if clust.name() in new_groups:
                    # The cluster has turned up in the new groups, so just update it
                    clust.reset_seq_ids(new_groups[clust.name()])
                    clust.sim_scores = clust.get_base_cluster().pull_scores_subgraph(new_groups[clust.name()])
                elif len(clust) > 1:
                    # The cluster didn't make it into new_groups, so it must have been emptied out
                    clust.reset_seq_ids([])
                    clust.sim_scores = clust.pull_scores_subgraph([])
                elif list(clust.seq_ids)[0] not in new_groups["orphans"]:
                    # An otherwise orphaned sequence was placed in a new cluster
                    clust.reset_seq_ids([])
                    clust.sim_scores = clust.pull_scores_subgraph([])
                else:
                    # Can only be a previously orphaned sequence. Remove it from new_group orphans to retain clust name.
                    new_groups["orphans"].remove(list(clust.seq_ids)[0])

            # Place any new orphans into new groups
            clusts_needing_update = []
            for seq_id in new_groups["orphans"]:
                sim_scores = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
                parent = None
                for clust in copy_clusters:
                    if seq_id in clust.seq_ids:
                        parent = [c for c in self.clusters if c.name() == clust.name()][0]
                        clusts_needing_update.append(parent.name())
                        break
                if parent is None:
                    raise AttributeError("Unable to find parent cluster")  # This should never happen

                # Placing the orphan with new cluster name
                clust = Cluster([seq_id], sim_scores, parent=parent)
                clust.set_name()
                self.clusters.append(clust)

            # If new clusters were made, need to update the originals too
            for clust_name in set(clusts_needing_update):
                clust = [clust for clust in self.clusters if clust.name() == clust_name][0]
                clust = Cluster(clust.seq_ids, clust.sim_scores, parent=clust)
                clust.set_name()
                self.clusters.append(clust)
                clust_names = [c.name() for c in self.clusters]
                del self.clusters[clust_names.index(clust_name)]

            # Delete any empty clusters
            del_list = []
            for indx, clust in enumerate(self.clusters):
                if not len(clust):
                    del_list.append(indx)
                    # Remove the associated group files if the group doesn't exist any longer
                    shutil.rmtree(join(self.outdir, "alignments", clust.name()), ignore_errors=True)
                    shutil.rmtree(join(self.outdir, "hmm", clust.name()), ignore_errors=True)
                    shutil.rmtree(join(self.outdir, "mcmcmc", clust.name()), ignore_errors=True)
                    shutil.rmtree(join(self.outdir, "sim_scores", "%s.scores" % clust.name()), ignore_errors=True)

            for indx in sorted(del_list, reverse=True):
                del self.clusters[indx]

            if breakout:
                break

        with LOCK:
            self.tmp_file.write(log_output.read())
        return


class FwdScoreCorrelations(object):
    def __init__(self, seqbuddy, outdir):
        self.seqbuddy = seqbuddy
        self.outdir = outdir
        self.tmp_dir = br.TempDir()
        self.tmp_dir.subfile("seqs.fa")

        self.hmm_fwd_scores_path = join(self.outdir, "hmm", "hmm_fwd_scores.csv")
        self.hmm_fwd_scores = self.create_hmm_fwd_score_df()

        self.rsquares_df_path = join(self.outdir, "hmm", "rsquares_matrix.csv")
        self.rsquare_vals_df = self.create_fwd_score_rsquared_matrix()

    def _mc_fwd_back_run(self, rec, args):
        hmm_scores_file = args[0]
        hmm_path = join(self.outdir, "hmm", "%s.hmm" % rec.id)

        fwdback_output = Popen("%s %s %s" % (HMM_FWD_BACK, hmm_path, join(self.tmp_dir.path, self.tmp_dir.subfiles[0])),
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

    def create_hmm_fwd_score_df(self, force=False):
        # Create HMM forward-score dataframe from initial input --> p(seq|hmm)
        if os.path.isfile(self.hmm_fwd_scores_path) and not force:
            return pd.read_csv(self.hmm_fwd_scores_path)

        self.seqbuddy.write(join(self.tmp_dir.path, self.tmp_dir.subfiles[0]), out_format="fasta")
        hmm_scores_file = br.TempFile()
        br.run_multicore_function(self.seqbuddy.records, self._mc_fwd_back_run, [hmm_scores_file.path],
                                  max_processes=CPUS, quiet=True)
        self.hmm_fwd_scores = pd.read_csv(hmm_scores_file.path, header=None)
        self.hmm_fwd_scores.columns = ["hmm_id", "rec_id", "fwd_raw"]
        self.hmm_fwd_scores.to_csv(self.hmm_fwd_scores_path, index=False)
        return self.hmm_fwd_scores

    @staticmethod
    def _mc_rsquare_vals(recs, args):
        rec1, sub_recs = recs
        hmm_fwd_scores, tmp_rsquares_file = args
        rsquare_vals_df = pd.DataFrame(columns=["seq1", "seq2", "r_square"])
        for rec2 in sub_recs:
            fwd1 = hmm_fwd_scores.loc[hmm_fwd_scores.rec_id == rec1.id].sort_values(by="hmm_id").fwd_raw
            fwd2 = hmm_fwd_scores.loc[hmm_fwd_scores.rec_id == rec2.id].sort_values(by="hmm_id").fwd_raw
            corr = scipy.stats.pearsonr(fwd1, fwd2)
            comparison = pd.DataFrame(data=[[rec1.id, rec2.id, corr[0]**2]],
                                      columns=["seq1", "seq2", "r_square"])
            rsquare_vals_df = rsquare_vals_df.append(comparison, ignore_index=True)
        rsquare_vals_df = rsquare_vals_df.append(pd.DataFrame(data=[[rec1.id, rec1.id, 1.0]],
                                                 columns=["seq1", "seq2", "r_square"]), ignore_index=True)

        rsquare_vals_df = rsquare_vals_df.to_csv(path_or_buf=None, header=None, index=False, index_label=False)
        with LOCK:
            with open(tmp_rsquares_file, "a") as ofile:
                ofile.write(rsquare_vals_df)

    def create_fwd_score_rsquared_matrix(self, force=False):
        # Calculate all-by-all matrix of correlation coefficients on fwd scores among all sequences
        if os.path.isfile(self.rsquares_df_path) and not force:
            return pd.read_csv(self.rsquares_df_path)

        sub_recs = list(self.seqbuddy.records[1:])
        multicore_data = [None for _ in range(len(self.seqbuddy))]
        tmp_rsquares_file = br.TempFile()
        for indx, rec1 in enumerate(self.seqbuddy.records):
            multicore_data[indx] = (rec1, sub_recs)
            sub_recs = sub_recs[1:]
        br.run_multicore_function(multicore_data, self._mc_rsquare_vals, [self.hmm_fwd_scores, tmp_rsquares_file.path],
                                  max_processes=CPUS)

        self.rsquare_vals_df = pd.read_csv(tmp_rsquares_file.path, header=None)
        self.rsquare_vals_df.columns = ["seq1", "seq2", "r_square"]
        self.rsquare_vals_df = self.rsquare_vals_df.sort_values(by=["seq1", "seq2"]).reset_index(drop=True)
        self.rsquare_vals_df.to_csv(self.rsquares_df_path, index=False)
        return self.rsquare_vals_df


def argparse_init():
    # Catch any extra parameters passed into -algn_p so they play nice with argparse
    if '--align_params' in sys.argv:
        sys.argv[sys.argv.index('--generate_alignment')] = '-algn_p'
    if '-algn_p' in sys.argv:
        ga_indx = sys.argv.index('-algn_p')
        sys.argv[ga_indx + 1] = " %s" % sys.argv[ga_indx + 1].rstrip()

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

  \033[1m\033[91mAll  input sequences  should  be  homologous! RD-MCL  is \033[4mNOT\033[24m  intended for
  identifying  orthogroups  from  whole  genome/transcriptome  data.\033[m

\033[1mUsage\033[m:
  rdmcl "/path/to/sequence_file" [-options]
''')

    parser.register('action', 'setup', _SetupAction)

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("sequences", help="Location of sequence file (most formats are fine)", action="store")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("-o", "--outdir", action="store", nargs="?", metavar="path",
                              help="Where should results be written? (default=rdmcl-%s-UID)" %
                                   time.strftime("%d-%m-%Y"))
    parser_flags.add_argument("-rs", "--r_seed", type=int, metavar="",
                              help="Specify a random seed for repeating a specific run")
    parser_flags.add_argument("-sql", "--sqlite_db", action="store", metavar="",
                              help="If RD-MCL has already created a SQLite database, reusing it can speed things up")
    parser_flags.add_argument("-psi", "--psipred_dir", action="store", metavar="",
                              help="If RD-MCL has already calculated PSI-Pred files, point to the directory")
    parser_flags.add_argument("-gn", "--group_name", action="store", default="group", metavar="",
                              help="Supply a name prefix for cluster names (default='group')")
    parser_flags.add_argument("-ts", "--taxa_sep", action="store", default="-", metavar="",
                              help="Specify the string that separates taxa ids from gene names (default='-')")
    parser_flags.add_argument("-ch", "--chains", default=MCMC_CHAINS, type=int, metavar="",
                              help="Specify how many MCMC chains to run (default=3)")
    parser_flags.add_argument("-wlk", "--walkers", default=3, type=int, metavar="",
                              help="Specify how many Metropolis-Hastings walkers are in each chain (default=3)")
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
    parser_flags.add_argument("-wdb", "--workdb", action="store", default="", metavar="",
                              help="Specify the directory that independent workers will be monitoring")
    parser_flags.add_argument("-algn_m", "--align_method", action="store", default="clustalo", metavar="",
                              help="Specify which alignment algorithm to use (supply full path if not in $PATH)")
    parser_flags.add_argument("-algn_p", "--align_params", action="store", default="", metavar="",
                              help="Supply alignment specific parameters")
    parser_flags.add_argument("-r", "--resume", action="store_true",
                              help="Try to pick up where a previous run left off (this breaks r_seed).")
    parser_flags.add_argument("-f", "--force", action="store_true",
                              help="Try to run no matter what.")
    parser_flags.add_argument("-q", "--quiet", action="store_true",
                              help="Suppress output during run (only final output is returned)")

    # Developer testing
    dev_flags = parser.add_argument_group(title="\033[1mDeveloper commands (Caution!)\033[m")
    dev_flags.add_argument("-spc", "--suppress_paralog_collapse", action="store_true",
                           help="Do not merge best hit paralogs")
    # dev_flags.add_argument("-sr", "--suppress_recursion", action="store_true",
    #                       help="Stop after a single round of MCL")
    # dev_flags.add_argument("-scc", "--suppress_clique_check", action="store_true",
    #                       help="Do not check for or break up cliques")
    dev_flags.add_argument("-ssf", "--suppress_singlet_folding", action="store_true",
                           help="Do not check for or merge singlets")
    # dev_flags.add_argument("-sit", "--suppress_iteration", action="store_true",
    #                       help="Only check for cliques and orphans once")
    dev_flags.add_argument("-trm", "--trimal", action="append", nargs="+", metavar="threshold",
                           help="Specify a list of trimal thresholds to apply (move from more strict to less)")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument("-setup", action="setup", dest=argparse.SUPPRESS, default=argparse.SUPPRESS)

    in_args = parser.parse_args()
    if not in_args.outdir:
        with open(in_args.sequences, "r") as ifile:
            uid = hlp.md5_hash(ifile.read() + str(in_args))[:6]
        in_args.outdir = join(os.getcwd(), "rdmcl-%s-%s" % (time.strftime("%d-%m-%Y"), uid))
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
    tmp_name = "".join([choice(br.string.ascii_letters + br.string.digits) for _ in range(10)])
    logger_obj = hlp.Logger(tmp_name)
    logging.info("*************************** Recursive Dynamic Markov Clustering ****************************")
    logging.warning(str(VERSION))
    logging.info("********************************************************************************************\n")

    if not os.path.isfile(join(SCRIPT_PATH, "config.ini")):
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

    if not shutil.which(in_args.align_method):
        logging.error("The alignment program '%s' is not detected "
                      "on your system." % in_args.align_method)
        sys.exit()

    # I'm pretty sure that SeqBuddy 1.2 is still fine, but Alb needs to be updated. Might as well just force the update.
    if int(re.search("([0-9]+)", Sb.VERSION.minor).group(1)) < 3:
        logging.error("\nERROR: Your version of SeqBuddy (%s.%s) is too old to work with RDMCL, please update it."
                      % (Sb.VERSION.major, Sb.VERSION.minor))
        sys.exit()
    logging.info("SeqBuddy version: %s.%s" % (Sb.VERSION.major, Sb.VERSION.minor))

    if int(re.search("([0-9]+)", Alb.VERSION.minor).group(1)) < 3:
        logging.error("Your version of AlignBuddy (%s.%s) is too old to work with RDMCL, please update it."
                      % (Alb.VERSION.major, Alb.VERSION.minor))
        sys.exit()
    logging.info("AlignBuddy version: %s.%s" % (Alb.VERSION.major, Sb.VERSION.minor))

    if not os.path.isfile(in_args.sequences):
        logging.error("Sequence file not found")
        sys.exit()

    sequences = Sb.SeqBuddy(in_args.sequences)
    sequences = Sb.clean_seq(sequences)  # Prevent any stray characters
    sequences = Sb.order_recs_by_len(sequences)  # Ensures all records are aligned on conserved domain during pileup
    tmp_alb = Alb.generate_msa(Sb.SeqBuddy(deepcopy(sequences.records[:2])), in_args.align_method, quiet=True)

    if tmp_alb.align_tool["tool"] == "MAFFT" and float(tmp_alb.align_tool["version"]) < 7.245:
        logging.error("\nERROR: Your version of MAFFT (%s) is too old to work with RDMCL, please update it."
                      % tmp_alb.align_tool["version"])
        sys.exit()

    global ALIGNMETHOD
    ALIGNMETHOD = in_args.align_method
    global ALIGNPARAMS
    ALIGNPARAMS = in_args.align_params

    if not ALIGNPARAMS:
        # If using clustal omega, mafft, or pagan, set the pileup flag
        tool = br.identify_msa_program(ALIGNMETHOD)
        if tool and tool["name"] in ["clustalo", "mafft"]:
            ALIGNPARAMS = "--pileup"
        elif tool and tool["name"] == "pagan":
            ALIGNPARAMS = "--pileup-alignment"

    logging.info("\nAlignment method: %s %s" % (tmp_alb.align_tool["tool"], tmp_alb.align_tool["version"]))
    logging.info("Alignment calls: $: %s %s" % (ALIGNMETHOD, ALIGNPARAMS))

    global TRIMAL
    if in_args.trimal:
        in_args.trimal = in_args.trimal[0]
        for indx, arg in enumerate(in_args.trimal):
            try:
                in_args.trimal[indx] = float(in_args.trimal[indx])
            except ValueError:
                pass
        TRIMAL = in_args.trimal
    logging.info("TrimAl values: %s" % TRIMAL)

    logging.warning("\nLaunching SQLite Daemons")

    sqlite_path = join(in_args.outdir, "sqlite_db.sqlite") if not in_args.sqlite_db \
        else os.path.abspath(in_args.sqlite_db)
    broker = hlp.SQLiteBroker(db_file=sqlite_path, lock_wait_time=in_args.lock_wait_time)
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()

    global WORKER_OUT
    global WORKER_DB
    global HEARTBEAT_DB

    if in_args.workdb:
        in_args.workdb = os.path.abspath(in_args.workdb)
        WORKER_OUT = join(in_args.workdb, ".worker_output")
        WORKER_DB = join(in_args.workdb, "work_db.sqlite")

        HEARTBEAT_DB = join(in_args.workdb, "heartbeat_db.sqlite")
        heartbeat = HeartBeat(HEARTBEAT_DB, MASTER_PULSE)
        heartbeat.start()
    else:
        heartbeat = HeartBeat("", "", dummy=True)

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
            heartbeat.end()
            sys.exit()
        else:
            logging.warning("Proceeding with large run.\n")

    if not check_sequences(sequences, in_args.taxa_sep):
        heartbeat.end()
        sys.exit()
    seq_ids_str = ", ".join(sorted([rec.id for rec in sequences.records]))
    seq_ids_hash = hlp.md5_hash(seq_ids_str)

    # Check if job has been run already
    if broker.query("SELECT (hash) FROM data_table WHERE hash=?", (seq_ids_hash,)) == seq_ids_hash:
        logging.warning("RESUME: This output directory was previous used for an identical RD-MCL run.\n"
                        "        All cached resources will be reused.")

    # Make sure all the necessary directories are present and emptied of old run files
    for _path in [join(in_args.outdir, x) for x in ["", "alignments", "mcmcmc", "sim_scores", "psi_pred", "hmm"]]:
        if not os.path.isdir(_path):
            logging.info("mkdir %s" % _path)
            os.makedirs(_path)
        # Delete old 'group' files/directories
        if not in_args.resume:
            root, dirs, files = next(os.walk(_path))
            for _file in files:
                if "group" in _file:
                    os.remove(join(root, _file))
            for _dir in dirs:
                if "group" in _dir:
                    shutil.rmtree(join(root, _dir))

    # Prepare log files into output directory
    if os.path.isfile(join(in_args.outdir, "rdmcl.log")):
        os.remove(join(in_args.outdir, "rdmcl.log"))

    logger_obj.move_log(join(in_args.outdir, "rdmcl.log"))

    if os.path.isfile(join(in_args.outdir, "placement.log")):
        os.remove(join(in_args.outdir, "placement.log"))

    if os.path.isfile(join(in_args.outdir, "cliques.log")):
        os.remove(join(in_args.outdir, "cliques.log"))

    with open(join(in_args.outdir, "paralog_cliques"), "w") as outfile:
        outfile.write("###########################################################\n"
                      "# If a named cluster contains reciprocal best hit cliques #\n"
                      "# among a group of paralogs, they are collapsed down to a #\n"
                      "# single representative. The collapses are itemized here. #\n"
                      "###########################################################\n\n")

    # Set CPU limits
    global CPUS
    CPUS = in_args.max_cpus

    # Write a copy of input sequences
    sequences.write(join(in_args.outdir, "input_seqs.fa"), out_format="fasta")

    # Remove medadata to reduce memory
    sequences = Sb.delete_metadata(sequences)

    # PSIPRED
    logging.warning("\n** PSI-Pred **")
    records_missing_ss_files = []
    records_with_ss_files = []
    if in_args.psipred_dir and not os.path.isdir(in_args.psipred_dir):
        logging.warning("PSI-Pred directory indicated below does not exits:\n%s \tUsing default location instead."
                        % in_args.psipred_dir)
        in_args.psipred_dir = ""

    in_args.psipred_dir = "{0}{1}psi_pred".format(in_args.outdir, os.sep) if not in_args.psipred_dir \
        else in_args.psipred_dir

    global PSIPREDDIR
    PSIPREDDIR = os.path.abspath(in_args.psipred_dir)

    for record in sequences.records:
        if os.path.isfile(join(in_args.psipred_dir, "%s.ss2" % record.id)):
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
        psi_pred_files.append((record.id, join(in_args.psipred_dir, "%s.ss2" % record.id)))

    psi_pred_files = OrderedDict(psi_pred_files)

    # Sequence HMMs
    logging.warning("\n** Creating Sequence-level HMMs **")
    br.run_multicore_function(sequences.records, mc_create_sequence_hmms, [in_args.outdir],
                              max_processes=CPUS, quiet=in_args.quiet)

    logging.warning("\n** Compile HMM forward-score correlation matrix **")
    fwd_scores_obj = FwdScoreCorrelations(sequences, in_args.outdir)

    # Initial alignment
    logging.warning("\n** All-by-all graph **")
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (in_args.open_penalty, in_args.ext_penalty))

    num_comparisons = ((len(sequences) ** 2) - len(sequences)) / 2
    logging.warning("Generating initial all-by-all similarity graph (%s comparisons)" % int(num_comparisons))
    logging.info(" written to: {0}{1}sim_scores{1}complete_all_by_all.scores".format(in_args.outdir, os.sep))
    scores_data, alignbuddy = retrieve_all_by_all_scores(sequences, psi_pred_files, broker)
    scores_data.to_csv(join(in_args.outdir, "sim_scores", "complete_all_by_all.scores"),
                       header=None, index=False, sep="\t")
    logging.info("\t-- finished in %s --\n" % TIMER.split())

    # Don't allow forward slashes or spaces in group names
    in_args.group_name = re.sub("[/ ]", "_", in_args.group_name)

    # First prepare the really raw first alignment in the database, without any collapsing.
    uncollapsed_group_0 = Cluster([rec.id for rec in sequences.records], scores_data,
                                  taxa_sep=in_args.taxa_sep, group_prefix=in_args.group_name, r_seed=in_args.r_seed)
    cluster2database(uncollapsed_group_0, broker, alignbuddy)

    # Then prepare the 'real' group_0 cluster
    if not in_args.suppress_paralog_collapse:
        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data, taxa_sep=in_args.taxa_sep,
                                  group_prefix=in_args.group_name, collapse=True, r_seed=in_args.r_seed,
                                  fwd_scores=fwd_scores_obj)
    else:
        group_0_cluster = uncollapsed_group_0

    if group_0_cluster.collapsed_genes:
        with open(join(in_args.outdir, "paralog_cliques"), "a") as outfile:
            outfile.write("# group_0\n")
            json.dump(group_0_cluster.collapsed_genes, outfile)
            outfile.write("\n\n")

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
                                       psi_pred_ss2=psi_pred_files, chains=MCMC_CHAINS, walkers=in_args.walkers,
                                       convergence=GELMAN_RUBIN, resume=in_args.resume)
    final_clusters = [cluster for cluster in final_clusters if cluster.subgroup_counter == 0]
    run_time.end()

    progress_dict = progress_tracker.read()
    logging.warning("Total MCL runs: %s" % progress_dict["mcl_runs"])
    logging.warning("\t-- finished in %s --" % TIMER.split())

    if not in_args.suppress_singlet_folding:
        logging.warning("\n** HMM-based sequence-to-cluster reassignment **")
        with open(join(in_args.outdir, "placement.log"), "a") as orphan_log_file:
            orphan_log_file.write("""\
    ######################################
    #  Initiating Sequence Reassignment  #
    ######################################

""")
        seq2clust_obj = Seqs2Clusters(final_clusters, 3, sequences, in_args.outdir)
        seq2clust_obj.place_seqs_in_clusts()
        with open(join(in_args.outdir, "placement.log"), "a") as orphan_log_file:
            orphan_log_file.write(seq2clust_obj.tmp_file.read())

        logging.warning("\t-- finished in %s --" % TIMER.split())
    # shutil.move("%s/%s" % (seq2clust_obj.tmp_dir.path, seq2clust_obj.tmp_dir.subdirs[0]), "hmms")
    # final_clusters = place_sequences_in_clusters(final_clusters, 3, sequences, in_args.outdir)
    '''
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
            with open(join(in_args.outdir, "cliques.log"), "a") as clique_log_file:
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
            with open(join(in_args.outdir, "placement.log"), "a") as orphan_log_file:
                orphan_log_file.write("""\
   #################################
   #  Initiating Orphan Placement  #
   #################################

""")
            orphans = Orphans(seqbuddy=sequences, clusters=final_clusters, sql_broker=broker,
                              psi_pred_ss2=psi_pred_files, outdir=in_args.outdir)
            orphans.place_orphans()
            final_clusters = orphans.clusters
            with open(join(in_args.outdir, "placement.log"), "a") as orphan_log_file:
                orphan_log_file.write(orphans.tmp_file.read())
            with open(join(in_args.outdir, "hmm", "pearsonr.csv"), "w") as pearsonr:
                orphans.rsquare_vals_df.to_csv(pearsonr)
            with open(join(in_args.outdir, "hmm", "fwd.csv"), "w") as fwd_file:
                orphans.hmm_fwd_scores.to_csv(fwd_file)

        if in_args.suppress_iteration:
            break

    if not in_args.suppress_clique_check or not in_args.suppress_singlet_folding:
        logging.warning("\t-- finished in %s --" % TIMER.split())
    '''

    logging.warning("\nPlacing any collapsed paralogs into their respective clusters")
    if group_0_cluster.collapsed_genes:
        with open(join(in_args.outdir, "paralog_cliques"), "a") as outfile:
            outfile.write("###########################################################\n"
                          "#  Paralogs were expanded back into the following groups  #\n"
                          "###########################################################\n\n")

    for clust in final_clusters:
        clust.parent = uncollapsed_group_0
        collapsed = False
        for seq_id in clust.seq_ids:
            if seq_id in group_0_cluster.collapsed_genes:
                clust.reset_seq_ids(group_0_cluster.collapsed_genes[seq_id] + list(clust.seq_ids))
                with open(join(in_args.outdir, "paralog_cliques"), "a") as outfile:
                    if not collapsed:
                        outfile.write("# %s\n" % clust.name())
                        collapsed = True
                    outfile.write('{"%s": %s}\n' % (seq_id, group_0_cluster.collapsed_genes[seq_id]))
        if collapsed:
            with open(join(in_args.outdir, "paralog_cliques"), "a") as outfile:
                outfile.write("\n")
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
            output += "\t".join(sorted(max_clust.seq_ids)) + '\n'
        del final_clusters[ind]

    logging.warning("\nTotal execution time: %s" % TIMER.total_elapsed())
    with open(join(in_args.outdir, "final_clusters.txt"), "w") as outfile:
        outfile.write(output)
        logging.warning("Final score: %s" % round(final_score, 4))
        if final_score < base_score:
            logging.warning("    This is lower than the initial base score after"
                            " the inclusion of collapsed paralogs")
        logging.warning("Clusters written to: %s%sfinal_clusters.txt" % (in_args.outdir, os.sep))

    heartbeat.end()
    broker.close()


def main():
    in_args = argparse_init()
    if 'setup' in in_args:
        setup()
    else:
        full_run(in_args)
    return


if __name__ == '__main__':
    main()
