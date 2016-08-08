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
Convert an all-by-all scores matrix into groups. MCMCMC-driven MCL is the basic work horse for clustering, and then
the groups are refined further to include singletons and doublets and/or be broken up into cliques, where appropriate.
"""

# Std library
import sys
import os
import re
import shutil
import logging
import json
from time import time
from copy import copy
from subprocess import Popen, PIPE, check_output
from multiprocessing import Lock
from random import random
from math import log
from hashlib import md5
from collections import OrderedDict

# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo

# My packages
import mcmcmc  # Note: This is in ../utilities and sym-linked to python3.5/site-packages
import MyFuncs
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb

git_commit = check_output(['git', 'rev-parse', '--short', 'HEAD']).decode().strip()
git_commit = " (git %s)" % git_commit if git_commit else ""
VERSION = "0.1.alpha%s" % git_commit
NOTICE = '''\
Public Domain Notice
--------------------
This is free software; see the LICENSE for further details.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov
'''


class Cluster(object):
    def __init__(self, seq_ids, sim_scores, out_dir=None, parent=None, clique=False):
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
        """
        self.taxa = OrderedDict()
        self.sim_scores = sim_scores
        self.parent = parent

        if not parent and not out_dir:
            raise AttributeError("Cannot create a root Cluster object without specifying an output directory.")
        self.out_dir = out_dir if out_dir else parent.out_dir

        self._subgroup_counter = 0
        self._clique_counter = 0
        self.clique = clique
        self.cliques = None
        self.cluster_score = None
        self.collapsed_genes = OrderedDict()  # If paralogs are reciprocal best hits, collapse them
        self._name = None
        self.similarity_graphs = "%s/sim_scores/all_graphs" % self.out_dir
        for next_seq_id in seq_ids:
            taxa = next_seq_id.split("-")[0]
            self.taxa.setdefault(taxa, [])
            self.taxa[taxa].append(next_seq_id)

        if clique and not parent:
            raise AttributeError("A clique cannot be declared without including its parental sequence_ids.")

        if parent:
            self.cluster_score_file = parent.cluster_score_file
            self.lock = parent.lock
            for indx, genes in parent.collapsed_genes.items():
                if indx in seq_ids:
                    self.collapsed_genes[indx] = genes
            self.seq_ids = seq_ids

        else:
            self._name = "group_0"
            # No reason to re-calculate scores, so record what's been done in a temp file to persist over subprocesses
            self.cluster_score_file = MyFuncs.TempFile()
            self.cluster_score_file.write("score,cluster")
            self.lock = Lock()

            # This next bit collapses all paralog reciprocal best-hit cliques so they don't gum up MCL
            breakout = False
            while not breakout:
                breakout = True
                for seq1_id in seq_ids:
                    seq1_taxa = seq1_id.split("-")[0]
                    paralog_best_hits = []
                    for hit in self._get_best_hits(seq1_id).itertuples():
                        seq2_id = hit.seq1 if hit.seq1 != seq1_id else hit.seq2
                        if seq2_id.split("-")[0] != seq1_taxa:
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
                        if paralog in self.collapsed_genes:
                            self.collapsed_genes[seq1_id] += self.collapsed_genes[paralog]
                            del self.collapsed_genes[paralog]

            self.seq_ids = seq_ids
            self.seq_id_hash = md5_hash("".join(sorted(self.seq_ids)))

    def _get_best_hits(self, gene):
        best_hits = self.sim_scores[(self.sim_scores.seq1 == gene) | (self.sim_scores.seq2 == gene)]
        if not best_hits.empty:
            best_hits = best_hits.loc[best_hits.score == max(best_hits.score)].values
            best_hits = pd.DataFrame(best_hits, columns=["seq1", "seq2", "score"])
        return best_hits

    def set_name(self):
        if self.name():
            pass
        elif not self.parent.name():
            raise ValueError("Parent of current cluster has not been named.\n%s" % self)
        elif self.clique:
            self._name = "%s_c%s" % (self.parent.name(), self.parent._clique_counter)
            self.parent._clique_counter += 1
        else:
            self._name = "%s_%s" % (self.parent.name(), self.parent._subgroup_counter)
            self.parent._subgroup_counter += 1
            if self.cliques and self.cliques[0]:
                for clique in self.cliques:
                    clique.set_name()
        return

    def compare(self, query):
        matches = set(self.seq_ids).intersection(query.seq_ids)
        weighted_match = (len(matches) * 2.) / (len(self) + query.len)
        print("name: %s, matches: %s, weighted_match: %s" % (self.name(), len(matches), weighted_match))
        return weighted_match

    ''' Delete this or refactor to do a new alignment for a proper subcluster (probably better to refactor).
    def sub_cluster(self, id_list):
        sim_scores = pd.DataFrame()
        for gene in id_list:
            if gene not in self.sequence_ids:
                raise NameError("Gene id '%s' not present in parental sequence_ids" % gene)
            seq1_scores = self.sim_scores[self.sim_scores.seq1 == gene]
            seq1_scores = seq1_scores[seq1_scores.seq2.isin(id_list)]
            seq2_scores = self.sim_scores[self.sim_scores.seq2 == gene]
            seq2_scores = seq2_scores[seq2_scores.seq1.isin(id_list)]
            scores = pd.concat([seq1_scores, seq2_scores])
            sim_scores = sim_scores.append(scores)
        sim_scores = sim_scores.drop_duplicates()
        subcluster = Cluster(id_list, sim_scores, _parent=self)
        return subcluster
    '''

    def name(self):
        return self._name

    def _recursive_best_hits(self, gene, global_best_hits, tested_ids):
        best_hits = self._get_best_hits(gene)
        global_best_hits = global_best_hits.append(best_hits, ignore_index=True)
        for _edge in best_hits.itertuples():
            if _edge.seq1 not in tested_ids:
                tested_ids.append(_edge.seq1)
                global_best_hits = self._recursive_best_hits(_edge.seq1, global_best_hits, tested_ids)
            if _edge.seq2 not in tested_ids:
                tested_ids.append(_edge.seq2)
                global_best_hits = self._recursive_best_hits(_edge.seq2, global_best_hits, tested_ids)
        return global_best_hits

    @staticmethod
    def perturb(scores):
        valve = MyFuncs.SafetyValve(global_reps=10)
        while scores.score.std() == 0:
            valve.step("Failed to perturb:\n%s" % scores)
            for indx, score in scores.score.iteritems():
                scores.set_value(indx, "score", random.gauss(score, (score * 0.0000001)))
        return scores

    def score(self):
        if self.cluster_score:
            return self.cluster_score
        # Confirm that a score for this cluster has not been calculated before
        prev_scores = pd.read_csv(self.cluster_score_file.path, index_col=False)
        seq_ids = md5_hash("".join(sorted(self.seq_ids)))
        if seq_ids in prev_scores.cluster.values:
            self.cluster_score = float(prev_scores.score[prev_scores.cluster == seq_ids])
            return self.cluster_score

        # Don't ignore the possibility of cliques, which will alter the score.
        # Note that cliques are assumed to be the smallest unit, so not containing any sub-cliques. Valid?
        if not self.cliques:
            best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
            # pull out all genes with replicate taxa and get best hits
            for taxa, genes in self.taxa.items():
                if len(genes) > 1:
                    for gene in genes:
                        best_hits = self._recursive_best_hits(gene, best_hits, [gene])

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
            self.cliques = []
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
                    clique = Cluster(clique, sim_scores=clique_scores, parent=self, clique=True)
                    self.cliques.append(clique)
            self.cliques = [None] if not self.cliques else self.cliques

        if self.cliques[0]:
            clique_list = [i for j in self.cliques for i in j.seq_ids]
            decliqued_cluster = []
            for gene in self.seq_ids:
                if gene not in clique_list:
                    decliqued_cluster.append(gene)
            decliqued_cluster = decliqued_cluster
        else:
            decliqued_cluster = self.seq_ids

        self.cluster_score = self.raw_score(decliqued_cluster)
        for clique in self.cliques:
            if not clique:
                break
            clique_ids = md5_hash("".join(sorted(clique.seq_ids)))
            if clique_ids in prev_scores.cluster.values:
                self.cluster_score += float(prev_scores.score[prev_scores.cluster == clique_ids])
            else:
                clique_score = self.raw_score(clique.seq_ids)
                self.cluster_score += clique_score
                with self.lock:
                    self.cluster_score_file.write("\n%s,%s" % (clique_score, clique_ids))
                    sim_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique.seq_ids)) &
                                                 (self.sim_scores.seq2.isin(clique.seq_ids))]
                    sim_scores.to_csv("%s/%s" % (self.similarity_graphs, clique_ids), index=False)
        with self.lock:
            self.cluster_score_file.write("\n%s,%s" % (self.cluster_score, seq_ids))
            self.sim_scores.to_csv("%s/%s" % (self.similarity_graphs, seq_ids), index=False)
        return self.cluster_score

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

    def raw_score(self, id_list):
        if len(id_list) == 1:
            return 0

        taxa = OrderedDict()
        for gene in id_list:
            taxon = gene.split("-")[0]
            taxa.setdefault(taxon, [])
            taxa[taxon].append(gene)

        # sequence_ids should never consist of a single taxa because reciprocal best hit paralogs have been removed
        # Update: Apparently not true. New paralog best hit cliques can form after the group is broken up some.
        # if len(taxa) == 1:
        #    raise ReferenceError("Only a single taxa found in sequence_ids...\n%s" % taxa)

        unique_scores = 1
        paralog_scores = -1
        for taxon, genes in taxa.items():
            if len(genes) == 1:
                unique_scores *= 2 * (1 + (len(self.taxa[taxon]) / len(self)))
            else:
                paralog_scores *= len(genes) * (len(genes) / len(self.taxa[taxon]))
        return unique_scores + paralog_scores

    def __len__(self):
        return len(self.seq_ids)

    def __str__(self):
        return str(self.seq_ids)


def md5_hash(in_str):
    in_str = str(in_str).encode()
    return md5(in_str).hexdigest()


def make_full_mat(subsmat):
    for key in copy(subsmat):
        try:
            # don't over-write the reverse keys if they are already initialized
            subsmat[(key[1], key[0])]
        except KeyError:
            subsmat[(key[1], key[0])] = subsmat[key]
    return subsmat


def bit_score(raw_score):
    # These values were empirically determined for BLOSUM62 by Altschul
    bit_k_value = 0.035
    bit_lambda = 0.252

    bits = ((bit_lambda * raw_score) - (log(bit_k_value))) / log(2)
    return bits


def psi_pred(seq_obj, args):
    out_dir = args[0]
    if os.path.isfile("%s/psi_pred/%s.ss2" % (out_dir, seq_obj.id)):
        return
    temp_dir = MyFuncs.TempDir()
    pwd = os.getcwd()
    psipred_dir = os.path.abspath("%s/psipred" % os.path.dirname(__file__))
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as ofile:
        ofile.write(seq_obj.format("fasta"))

    command = '''\
psiblast -db {0}/blastdb/pannexins -query sequence.fa -inclusion_ethresh 0.001 -out_pssm {1}/{2}.chk \
-num_iterations 3 -num_alignments 0 >& {1}/{2}.blast;
{0}/bin/chkparse {1}/{2}.chk > {1}/{2}.mtx;
{0}/bin/psipred {1}/{2}.mtx {0}/data/weights.dat {0}/data/weights.dat2 {0}/data/weights.dat3 > {1}/{2}.ss;
{0}/bin/psipass2 {0}/data/weights_p2.dat 1 1.0 1.0 {1}/{2}.ss2 {1}/{2}.ss > {1}/{2}.horiz;
'''.format(psipred_dir, temp_dir.path, seq_obj.id)

    Popen(command, shell=True).wait()
    os.chdir(pwd)
    shutil.move("%s/%s.ss2" % (temp_dir.path, seq_obj.id), "%s/psi_pred/" % out_dir)
    return


def mcmcmc_mcl(args, params):
    inflation, gq = args
    external_tmp_dir, min_score, seqbuddy, parent_cluster = params
    mcl_tmp_dir = MyFuncs.TempDir()

    mcl_output = Popen("mcl %s/input.csv --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
                       (external_tmp_dir, gq, inflation, mcl_tmp_dir.path), shell=True, stderr=PIPE).communicate()

    with lock:
        with open("%s/.progress" % in_args.outdir, "r") as ifile:
            finished_clusters, mcl_runs = ifile.read().split("\n")
        with open("%s/.progress" % in_args.outdir, "w") as ofile:
            ofile.write("%s\n%s" % (finished_clusters, int(mcl_runs) + 1))

    mcl_output = mcl_output[1].decode()
    if re.search("\[mclvInflate\] warning", mcl_output) and min_score:
        return min_score

    clusters = parse_mcl_clusters("%s/output.groups" % mcl_tmp_dir.path)
    score = 0
    prev_scores = pd.read_csv(parent_cluster.cluster_score_file.path, index_col=False)

    for indx, cluster in enumerate(clusters):
        cluster_ids = md5_hash("".join(sorted(cluster)))
        if cluster_ids in prev_scores.cluster.values:
            with parent_cluster.lock:
                sim_scores = pd.read_csv("%s/%s" % (parent_cluster.similarity_graphs, cluster_ids), index_col=False)
                scores = prev_scores.score[prev_scores.cluster == cluster_ids]
                if len(scores) > 1:
                    prev_scores = prev_scores.drop_duplicates()
                    prev_scores.to_csv(parent_cluster.cluster_score_file.path, index=False)
            score += float()
        else:
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in cluster]))
            alb_obj = generate_msa(sb_copy)
            sim_scores = create_all_by_all_scores(alb_obj, quiet=True)

        cluster = Cluster(cluster, sim_scores, parent=parent_cluster)
        clusters[indx] = cluster
        score += cluster.score()

    with lock:
        with open("%s/max.txt" % external_tmp_dir, "r") as ifile:
            current_max = float(ifile.read())
        if score > current_max:
            write_mcl_clusters(clusters, "%s/best_group" % external_tmp_dir)
            with open("%s/max.txt" % external_tmp_dir, "w") as ofile:
                ofile.write(str(score))
    return score


def orthogroup_caller(master_cluster, cluster_list, seqbuddy, steps=1000, quiet=True):
    """
    Run MCMCMC on MCL to find the best orthogroups
    :param master_cluster: The group to be subdivided
    :type master_cluster: Cluster
    :param cluster_list: When a sequence_ids is finalized after recursion, it is appended to this list
    :param seqbuddy: The sequences that are included in the master sequence_ids
    :param steps: How many MCMCMC iterations to run TODO: calculate this on the fly
    :param quiet: Suppress StdErr
    :return: list of sequence_ids objects
    """
    def save_cluster():
        if not os.path.isdir("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name())):
            temp_dir.save("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name()))
        alignment = generate_msa(seqbuddy)
        alignment.write("%s/alignments/%s.aln" % (in_args.outdir, master_cluster.name()))
        master_cluster.sim_scores.to_csv("%s/sim_scores/%s.scores" % (in_args.outdir, master_cluster.name()),
                                         header=None, index=False, sep="\t")
        return

    temp_dir = MyFuncs.TempDir()
    master_cluster.sim_scores.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.1, 20)
    gq_var = mcmcmc.Variable("gq", min(master_cluster.sim_scores.score), max(master_cluster.sim_scores.score))

    try:
        with open("%s/max.txt" % temp_dir.path, "w") as ofile:
            ofile.write("-1000000000")
        with lock:
            with open("%s/.progress" % in_args.outdir, "r") as ifile:
                finished_clusters, mcl_runs = ifile.read().split("\n")

            with open("%s/.progress" % in_args.outdir, "w") as ofile:
                ofile.write("%s\n%s" % (len(cluster_list), mcl_runs))

        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=["%s" % temp_dir.path, False, seqbuddy, master_cluster], quiet=quiet,
                                       outfile="%s/mcmcmc_out.csv" % temp_dir.path)

    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        save_cluster()
        return cluster_list

    # Set a 'worst score' that is reasonable for the data set
    worst_score = 10000000  # arbitrarily large number to start
    for chain in mcmcmc_factory.chains:
        worst_score = chain.raw_min if chain.raw_min < worst_score else worst_score

    mcmcmc_factory.reset_params(["%s" % temp_dir.path, worst_score, seqbuddy, master_cluster])
    mcmcmc_factory.run()
    mcmcmc_output = pd.read_csv("%s/mcmcmc_out.csv" % temp_dir.path, "\t")

    best_score = max(mcmcmc_output["result"])
    if best_score <= master_cluster.score():
        if master_cluster.name() != "group_0":
            master_cluster.set_name()
            cluster_list.append(master_cluster)
            save_cluster()
        return cluster_list

    mcl_clusters = parse_mcl_clusters("%s/best_group" % temp_dir.path)
    for sub_cluster in mcl_clusters:
        cluster_ids = md5_hash("".join(sorted(sub_cluster)))
        sim_scores = pd.read_csv("%s/%s" % (master_cluster.similarity_graphs, cluster_ids), index_col=False)
        sub_cluster = Cluster(sub_cluster, sim_scores=sim_scores, parent=master_cluster)
        if len(sub_cluster) in [1, 2]:
            sub_cluster.set_name()
            cluster_list.append(sub_cluster)
            continue

        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in sub_cluster.seq_ids])

        # Recursion...
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, steps=steps, quiet=quiet,)

    save_cluster()
    return cluster_list


def progress():
    try:
        with open("%s/.progress" % in_args.outdir, "r") as ifile:
            finished_clusters, mcl_runs = ifile.read().split("\n")
        return "MCL runs processed: %s. Clusters finished: %s. Run time: " % (mcl_runs, finished_clusters)
    except ValueError:  # In case the file is being written while trying to read
        return "MCL runs processed: ?. Clusters finished: ?. Run time: "


def parse_mcl_clusters(path):
    with open(path, "r") as ifile:
        clusters = ifile.read()
    clusters = clusters.strip().split("\n")
    clusters = [cluster.strip().split("\t") for cluster in clusters]
    return clusters


def write_mcl_clusters(clusters, path):
    clusters_string = ""
    for cluster in clusters:
        clusters_string += "%s\n" % "\t".join(cluster.seq_ids)
    with open(path, "w") as ofile:
        ofile.write(clusters_string.strip())
    return


def place_orphans(clusters, scores):
    """
    If a cluster has only one or two sequences in it, it may make sense to fold that cluster into one of the larger
    clusters. To check if this is reasonable, compare the similarity scores between the orphan sequences and the
    sequences in each larger cluster. Perform a Tukey HSD test among all of these comparisons and see if one cluster
    has significantly higher similarity scores than ALL of the other clusters.
    NOTE: Placing an orphan into a cluster obviously changes the content of that larger cluster, but this change is
    ignored. Perhaps this function could be run repeatedly until there are no more changes, but for now it's only
    called once.
    :param clusters: The entire complement of clusters returned by RD-MCL (list of Cluster objects)
    :param scores: The initial all-by-all similarity graph used for RD-MCL
    :return: Updated clusters list
    """
    small_clusters = OrderedDict()
    large_clusters = OrderedDict()
    for cluster in clusters:
        if len(cluster.seq_ids) > 2:
            large_clusters[cluster.name()] = cluster
        else:
            small_clusters[cluster.name()] = cluster

    if not small_clusters:
        logging.warning("No orphaned sequences present")
        return clusters
    elif not large_clusters:
        logging.warning("All clusters have only 1 or 2 sequences present, suggesting there are issues with your data")
        logging.info(" Are there a large number of paralogs and only a small number of taxa present?")
        return clusters
    elif len(large_clusters) == 1:
        logging.warning("Only one orthogroup with >2 sequences present, suggesting there are issues with your data")
        logging.info(" Are there a large number of taxa and only a small number of orthologs present?")
        return clusters
    else:
        num_orphans = sum([len(cluster) for group_name, cluster in small_clusters.items()])
        logging.warning("%s orphaned sequences present" % num_orphans)

    # Create lists of similarity scores for each orphan cluster against each larger cluster
    def get_score(seq1, seq2):
        sim_score = scores.loc[:][scores.seq1 == seq1]
        sim_score = sim_score.loc[:][sim_score.seq2 == seq2]

        if sim_score.empty:
            sim_score = scores.loc[:][scores.seq1 == seq2]
            sim_score = sim_score.loc[:][sim_score.seq2 == seq1]

        if sim_score.empty:
            sim_score = 0.
        else:
            sim_score = float(sim_score.score)
        return sim_score

    large_clusters_ordered_dict = OrderedDict([(group_name, []) for group_name in large_clusters])
    small_to_large_map = [(group_name, large_clusters_ordered_dict.copy()) for group_name in small_clusters]
    small_to_large_map = OrderedDict(small_to_large_map)
    for sgroup, sclust in small_clusters.items():           # Get small cluster
        for sgene in sclust.seq_ids:                        # Get gene id from small cluster
            for lgroup, lclust in large_clusters.items():   # Get large cluster
                for lgene in lclust.seq_ids:                # Get gene id from large cluster
                    small_to_large_map[sgroup][lgroup].append(get_score(sgene, lgene))  # Update appropriate list

    for sgroup, l_clusts in small_to_large_map.items():
        # Convert data into list of numpy arrays that sm.stats can read, also get average scores for each sequence_ids
        data = [np.array(x) for x in l_clusts.values()]

        # Only need to test the cluster with the highest average similarity scores, so find which cluster that is.
        averages = pd.Series()
        for j, group in enumerate(data):
            key = list(l_clusts.keys())[j]
            averages = averages.append(pd.Series(np.mean(group), index=[key]))
            data[j] = pd.DataFrame(group, columns=['observations'])
            data[j]['grouplabel'] = key
        max_ave = averages.argmax()
        df = pd.DataFrame()
        for group in data:
            df = df.append(group)

        # Run pairwise Tukey HSD and parse the results
        result = sm.stats.multicomp.pairwise_tukeyhsd(df.observations, df.grouplabel)
        for line in str(result).split("\n")[4:-1]:
            line = re.sub("^ *", "", line.strip())
            line = re.sub(" +", "\t", line)
            line = line.split("\t")
            if max_ave in line and 'False' in line:
                break  # Insufficient support to group the gene with max_ave group

        else:  # If the 'break' command is not encountered, the gene can be grouped with the max_ave cluster
            large_clusters[max_ave].seq_ids += small_clusters[sgroup].seq_ids
            logging.info("\t%s added to %s" % (" and ".join(small_clusters[sgroup].seq_ids), max_ave))
            del small_clusters[sgroup]

    clusters = [cluster for key, cluster in large_clusters.items()]
    clusters += [cluster for key, cluster in small_clusters.items()]
    fostered_orphans = num_orphans - sum([len(cluster) for key, cluster in small_clusters.items()])
    if not fostered_orphans:
        logging.warning("No orphaned sequences were placed in orthogroups.")
    elif fostered_orphans == 1:
        logging.warning("1 orphaned sequence was placed in an orthogroups.")
    else:
        logging.warning("%s orphaned sequences have found new homes in orthogroups!" % fostered_orphans)
        logging.warning("\tNote that all of the saved alignment, mcmcmc, and sim_score\n"
                        "\tfiles will only include the original sequences.")
    return clusters
# NOTE: There used to be a support function. Check the workscripts GitHub history


def score_sequences(seq_pair, args):
    # Calculate the best possible scores, and divide by the observed scores
    id1, id2 = seq_pair
    alb_obj, psi_pred_files, putput_file = args
    id_regex = "^%s$|^%s$" % (id1, id2)
    alb_copy = Alb.make_copy(alb_obj)
    Alb.pull_records(alb_copy, id_regex)
    observed_score = 0
    seq1_best = 0
    seq2_best = 0
    seq1, seq2 = alb_copy.records()
    prev_aa1 = "-"
    prev_aa2 = "-"

    for aa_pos in range(alb_copy.lengths()[0]):
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

    subs_mat_score = ((observed_score / seq1_best) + (observed_score / seq1_best)) / 2

    # PSI PRED comparison
    num_gaps = 0
    ss_score = 0
    for row1 in psi_pred_files[id1].itertuples():
        if (psi_pred_files[id2]["indx"] == row1.indx).any():
            row2 = psi_pred_files[id2].loc[psi_pred_files[id2]["indx"] == row1.indx]
            row_score = 0
            row_score += 1 - abs(float(row1.coil_prob) - float(row2.coil_prob))
            row_score += 1 - abs(float(row1.helix_prob) - float(row2.helix_prob))
            row_score += 1 - abs(float(row1.sheet_prob) - float(row2.sheet_prob))
            ss_score += row_score / 3
        else:
            num_gaps += 1

    align_len = len(psi_pred_files[id2]) + num_gaps
    ss_score /= align_len
    final_score = (ss_score * 0.3) + (subs_mat_score * 0.7)
    with lock:
        with open(putput_file, "a") as ofile:
            ofile.write("\n%s,%s,%s" % (id1, id2, final_score))
    return


def generate_msa(seqbuddy):
    seq_ids = sorted([rec.id for rec in seqbuddy.records])
    seq_id_hash = md5_hash("".join(seq_ids))
    align_file = "%s/alignments/all_alignments/%s" % (in_args.outdir, seq_id_hash)
    if os.path.isfile(align_file):
        alignment = Alb.AlignBuddy(align_file)
    else:
        if len(seqbuddy) == 1:
            alignment = Alb.AlignBuddy(str(seqbuddy))
        else:
            alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft", params="--globalpair --thread -1", quiet=True)
        alignment.write(align_file)
    return alignment


def create_all_by_all_scores(alignment, quiet=False):
    """
    Generate a multiple sequence alignment and pull out all-by-all similarity graph
    :param alignment: AlignBuddy object
    :param quiet: Supress multicore output
    :return:
    """
    if len(alignment.records()) == 1:
        sim_scores = pd.DataFrame(data=None, columns=["seq1", "seq2", "score"])
        return sim_scores

    # Only calculate if not previously calculated
    seq_ids = sorted([rec.id for rec in alignment.records_iter()])
    sim_scores_file = "%s/sim_scores/all_graphs/%s" % (in_args.outdir, md5_hash("".join(seq_ids)))
    if os.path.isfile(sim_scores_file):
        sim_scores = pd.read_csv(sim_scores_file, index_col=False)
        sim_scores.columns = ["seq1", "seq2", "score"]
        return sim_scores

    # Don't want to modify the alignbuddy object in place
    alignment = Alb.make_copy(alignment)

    # Need to specify what columns the PsiPred files map to now that there are gaps.
    psi_pred_files = OrderedDict()
    for rec in alignment.records_iter():
        ss_file = pd.read_csv("%s/psi_pred/%s.ss2" % (in_args.outdir, rec.id), comment="#",
                              header=None, delim_whitespace=True)
        ss_file.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
        ss_counter = 0
        for indx, residue in enumerate(rec.seq):
            if residue != "-":
                ss_file.set_value(ss_counter, "indx", indx)
                ss_counter += 1
        psi_pred_files[rec.id] = ss_file

    # Scored seem to be improved by removing gaps. Need to test this explicitly for the paper though
    alignment = Alb.trimal(alignment, "gappyout")

    # Re-update PsiPred files, now that some columns are removed
    for rec in alignment.records_iter():
        new_psi_pred = []
        for row in psi_pred_files[rec.id].itertuples():
            if alignment.alignments[0].position_map[int(row[1])][1]:
                new_psi_pred.append(list(row)[1:])
        psi_pred_files[rec.id] = pd.DataFrame(new_psi_pred, columns=["indx", "aa", "ss", "coil_prob",
                                                                     "helix_prob", "sheet_prob"])
    ids1 = [rec.id for rec in alignment.records_iter()]
    ids2 = [rec.id for rec in alignment.records_iter()]
    all_by_all = []
    for rec1 in ids1:
        del ids2[ids2.index(rec1)]
        for rec2 in ids2:
            all_by_all.append((rec1, rec2))

    ofile = MyFuncs.TempFile()
    ofile.write("seq1,seq2,score")
    printer.clear()
    MyFuncs.run_multicore_function(all_by_all, score_sequences, [alignment, psi_pred_files, ofile.path], quiet=quiet)
    sim_scores = pd.read_csv(ofile.path, index_col=False)
    return sim_scores


class Logger(object):
    def __init__(self, location=None):
        if not location:
            tmpfile = MyFuncs.TempFile()
            self.location = "%s/rdmcl.log" % tmpfile.path
        else:
            self.location = location

        # Set up logging. Use 'info' to write to file only, anything higher will go to both terminal and file.
        logging.basicConfig(filename=location, level=logging.INFO, format="")
        self.logger = logging.getLogger()
        self.console = logging.StreamHandler()
        self.console.setLevel(logging.WARNING)
        self.logger.addHandler(self.console)

    def move_log(self, location):
        shutil.move(self.location, location)
        logging.basicConfig(filename=location, level=logging.INFO, format="")
        self.location = location
        return


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(prog="orthogroup_caller", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("sequences", help="Location of sequence file", action="store")
    parser.add_argument("outdir", action="store", default="%s/rd-mcd" % os.getcwd(),
                        help="Where should results be written?")
    parser.add_argument("-sz", "--sample_size", type=float, default=0.632,
                        help="Proportion of total population to use in each jackknife replicate")
    parser.add_argument("-mcs", "--mcmcmc_steps", default=1000, type=int,
                        help="Specify how deeply to sample MCL parameters")
    parser.add_argument("-sr", "--supress_recursion", action="store_true",
                        help="Stop after a single round of MCL. For testing.")
    parser.add_argument("-scc", "--supress_clique_check", action="store_true",
                        help="Do not check for or break up cliques. For testing.")
    parser.add_argument("-ssf", "--supress_singlet_folding", action="store_true",
                        help="Do not check for or merge singlets. For testing.")
    parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                        type=float, default=-5)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=0)
    parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                        help="Don't apply the gappyout algorithm to MSAs before scoring")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Overwrite previous run")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress all output during run (only final output is returned)")

    in_args = parser.parse_args()

    logger_obj = Logger("rdmcl.log")
    logging.info("*************************** Recursive Dynamic Markov Clustering ****************************")
    logging.warning("RD-MCL version %s\n\n%s" % (VERSION, NOTICE))
    logging.info("********************************************************************************************\n")
    logging.info("Working directory: %s" % os.getcwd())
    logging.info("Function call: %s" % " ".join(sys.argv))

    sequences = Sb.SeqBuddy(in_args.sequences)
    seq_ids_hash = md5_hash("".join(sorted([rec.id for rec in sequences.records])))
    previous_run = False

    if os.path.exists("%s/.rdmcl" % in_args.outdir):
        with open("%s/.rdmcl" % in_args.outdir, "r") as hash_file:
            previous_run = hash_file.read()

    if os.path.exists(in_args.outdir):
        if previous_run == seq_ids_hash:
            logging.warning("RESUME: This output directory was previous used for an identical RD-MCL run; any\n"
                            "        cached resources will be reused.")
        else:
            check = MyFuncs.ask("Output directory already exists, continue [y]/n?") if not in_args.force else True
            if not check:
                logging.warning("Program aborted by user.")
                sys.exit()

    # Make sure all the necessary directories are present and emptied of old run files
    for outdir in ["%s%s" % (in_args.outdir, x) for x in ["", "/alignments", "/mcmcmc", "/sim_scores", "/psi_pred"]]:
        if not os.path.isdir(outdir):
            logging.info("mkdir %s" % outdir)
            os.makedirs(outdir)
        elif "psi_pred" not in outdir:  # Delete old files
            root, dirs, files = next(os.walk(outdir))
            for file in files:
                os.remove("%s/%s" % (root, file))

    with open("%s/.rdmcl" % in_args.outdir, "w") as hash_file:
        hash_file.write(seq_ids_hash)

    if not os.path.isdir("%s/sim_scores/all_graphs" % in_args.outdir):
        os.makedirs("%s/sim_scores/all_graphs" % in_args.outdir)

    if not os.path.isdir("%s/alignments/all_alignments" % in_args.outdir):
        os.makedirs("%s/alignments/all_alignments" % in_args.outdir)

    root, dirs, files = next(os.walk("%s/mcmcmc" % in_args.outdir))
    for _dir in dirs:
        shutil.rmtree("%s/%s" % (root, _dir))

    # Move log file into output directory
    logger_obj.move_log("%s/rdmcl.log" % in_args.outdir)

    clique_check = True if not in_args.supress_clique_check else False
    recursion_check = True if not in_args.supress_recursion else False
    best = None
    best_clusters = None
    lock = Lock()
    printer = MyFuncs.DynamicPrint(quiet=in_args.quiet)
    timer = MyFuncs.Timer()

    BLOSUM62 = make_full_mat(SeqMat(MatrixInfo.blosum62))

    ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
                   "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
    for aa in ambiguous_X:
        pair = sorted((aa, "X"))
        pair = tuple(pair)
        BLOSUM62[pair] = ambiguous_X[aa]

    # PSIPRED
    logging.warning("\n** PSI-Pred **")
    records_missing_ss_files = []
    records_with_ss_files = []
    for record in sequences.records:
        if os.path.isfile("%s/psi_pred/%s.ss2" % (in_args.outdir, record.id)):
            records_with_ss_files.append(record.id)
        else:
            records_missing_ss_files.append(record)
    if records_missing_ss_files and len(records_missing_ss_files) != len(sequences):
        logging.info("RESUME: PSI-Pred .ss2 files found for %s sequences:" % len(records_with_ss_files))

    if records_missing_ss_files:
        logging.warning("Executing PSI-Pred on %s sequences" % len(records_missing_ss_files))
        MyFuncs.run_multicore_function(records_missing_ss_files, psi_pred, [in_args.outdir])
        logging.info("\t-- finished in %s --" % timer.split())
        logging.info("\tfiles saved to %s\n" % "%s/psi_pred/" % in_args.outdir)
    else:
        logging.warning("RESUME: All PSI-Pred .ss2 files found in %s/psi_pred/" % in_args.outdir)

    # Initial alignment
    logging.warning("\n** All-by-all graph **")
    gap_open = in_args.open_penalty
    gap_extend = in_args.extend_penalty
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (gap_open, gap_extend))

    if os.path.isfile("%s/alignments/all_alignments/%s" % (in_args.outdir, seq_ids_hash)):
        logging.warning("RESUME: Initial multiple sequence alignment found")
        alignbuddy = Alb.AlignBuddy("%s/alignments/all_alignments/%s" % (in_args.outdir, seq_ids_hash))
    else:
        logging.warning("Generating initial multiple sequence alignment with MAFFT")
        alignbuddy = generate_msa(sequences)
        alignbuddy.write("%s/alignments/group_0.aln" % in_args.outdir)
        alignbuddy.write("%s/alignments/all_alignments/%s" % (in_args.outdir, seq_ids_hash))
        logging.info("\t-- finished in %s --" % timer.split())

    if os.path.isfile("%s/sim_scores/all_graphs/%s" % (in_args.outdir, seq_ids_hash)):
        logging.warning("RESUME: Initial all-by-all similarity graph found")
        scores_data = pd.read_csv("%s/sim_scores/all_graphs/%s" % (in_args.outdir, seq_ids_hash), index_col=False)
        scores_data.columns = ["seq1", "seq2", "score"]
        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data, out_dir=in_args.outdir)

    else:
        logging.warning("Generating initial all-by-all similarity graph")
        scores_data = create_all_by_all_scores(alignbuddy)
        scores_data.to_csv("%s/sim_scores/group_0.scores" % in_args.outdir, header=None, index=False, sep="\t")
        scores_data.to_csv("%s/sim_scores/all_graphs/%s" % (in_args.outdir, seq_ids_hash), header=None, index=False)
        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data, out_dir=in_args.outdir)
        logging.info("\t-- finished in %s --" % timer.split())

    # Base cluster score
    base_score = group_0_cluster.score()
    logging.warning("Base cluster score: %s" % round(base_score, 4))
    if group_0_cluster.collapsed_genes:
        logging.warning("Reciprocal best-hit cliques of paralogs have been identified in the input sequences.")
        logging.info(" A representative sequence has been selected from each clique, and the remaining")
        logging.info(" sequences will be placed in the final clusters at the end of the run.")
        with open("%s/paralog_cliques.json" % in_args.outdir, "w") as outfile:
            json.dump(group_0_cluster.collapsed_genes, outfile)
            logging.warning(" Cliques written to: %s/paralog_cliques.json" % in_args.outdir)

    # taxa_count = [x.split("-")[0] for x in group_0_cluster.seq_ids]
    # taxa_count = pd.Series(taxa_count)
    # taxa_count = taxa_count.value_counts()

    # Ortholog caller
    logging.warning("\n** Recursive MCL **")
    final_clusters = []
    with open("%s/.progress" % in_args.outdir, "w") as progress_file:
        progress_file.write("0\n0")
    run_time = MyFuncs.RunTime(prefix=progress, sleep=0.3)
    run_time.start()
    final_clusters = orthogroup_caller(group_0_cluster, final_clusters, seqbuddy=sequences,
                                       steps=in_args.mcmcmc_steps, quiet=True)
    run_time.end()
    with open("%s/.progress" % in_args.outdir, "r") as progress_file:
        progress_file = progress_file.read().split("\n")
    logging.info("Total MCL runs: %s" % progress_file[1])
    logging.info("\t-- finished in %s --" % timer.split())

    # Fold singletons and doublets back into groups. This can't be 'resumed', because it changes the clusters
    if not in_args.supress_singlet_folding:
        logging.warning("\n** Folding orphan sequences into clusters **")
        final_clusters = place_orphans(final_clusters, group_0_cluster.sim_scores)
        logging.warning("\t-- finished in %s --" % timer.split())

    # Format the clusters and output to file
    logging.warning("\n** Final formatting **")
    if group_0_cluster.collapsed_genes:
        logging.warning("Placing collapsed paralogs into their respective clusters")
        for clust in final_clusters:
            for gene_id, paralogs in group_0_cluster.collapsed_genes.items():
                if gene_id in clust.seq_ids:
                    clust.seq_ids += paralogs

    logging.warning("Preparing final_clusters.txt")
    output = ""
    while len(final_clusters) > 0:
        _max = (0, 0)
        for ind, clust in enumerate(final_clusters):
            if len(clust.seq_ids) > _max[1]:
                _max = (ind, len(clust.seq_ids))

        ind, max_clust = _max[0], final_clusters[_max[0]]
        del final_clusters[ind]
        output += "group_%s\t%s\t" % (max_clust.name(), round(max_clust.score(), 4))
        for seq_id in max_clust.seq_ids:
            output += "%s\t" % seq_id
        output = "%s\n" % output.strip()

    logging.warning("\t-- finished in %s --" % timer.split())

    logging.warning("\nTotal execution time: %s" % timer.total_elapsed())
    with open("%s/final_clusters.txt" % in_args.outdir, "w") as outfile:
        outfile.write(output)
        logging.warning("Final clusters written to: %s/final_clusters.txt" % in_args.outdir)
