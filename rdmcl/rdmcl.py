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
from time import time
from copy import copy
from subprocess import Popen, PIPE, check_output
from multiprocessing import Lock
from random import random
from math import log
from hashlib import md5

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
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov
'''


class Cluster(object):
    def __init__(self, seq_ids, sim_scores, _parent=None, clique=False):
        """
        - Note that reciprocal best hits between paralogs are collapsed when instantiating group_0, so
          no problem strongly penalizing all paralogs in the scoring algorithm

        :param seq_ids: Sequence IDs
        :type seq_ids: list
        :param sim_scores: All-by-all similarity matrix for everything in the seq_ids (and parental clusters)
        :type sim_scores: pandas.DataFrame
        :param _parent: Parental seq_ids
        :type _parent: Cluster
        :param clique: Specify whether cluster is a clique or not
        :type clique: bool
        """
        self.taxa = {}
        self.sim_scores = sim_scores
        self.parent = _parent
        self._subgroup_counter = 0
        self._clique_counter = 0
        self.clique = clique
        self.cliques = None
        self._score = None
        self.collapsed_genes = {}  # If paralogs are reciprocal best hits, collapse them
        self._name = None

        if clique and not _parent:
            raise AttributeError("A clique cannot be declared without including its parental seq_ids.")

        if _parent:
            self.cluster_score_file = _parent.cluster_score_file
            self.similarity_graphs = _parent.similarity_graphs
            self.lock = _parent.lock
            for indx, genes in _parent.collapsed_genes.items():
                if indx in seq_ids:
                    self.collapsed_genes[indx] = genes
            self.seq_ids = seq_ids
            for gene in seq_ids:
                taxa = gene.split("-")[0]
                self.taxa.setdefault(taxa, [])
                self.taxa[taxa].append(gene)
        else:
            self._name = "group_0"
            # No reason to re-calculate scores, so record what's been done in a temp file to persist over subprocesses
            self.cluster_score_file = MyFuncs.TempFile()
            self.cluster_score_file.write("score,cluster")
            self.similarity_graphs = MyFuncs.TempDir()
            self.lock = Lock()
            collapse_check_list = []
            collapsed_cluster = []
            for gene in seq_ids:
                if gene in collapse_check_list:
                    continue

                taxa = gene.split("-")[0]
                self.taxa.setdefault(taxa, [])
                self.taxa[taxa].append(gene)
                collapsed_cluster.append(gene)
                valve = MyFuncs.SafetyValve()  # ToDo: remove this when sure it is unnecessary
                breakout = False
                while not breakout:  # This is the primary paralog collapsing logic. Only do it once for parental group
                    valve.step()
                    breakout = True
                    for edge in self._get_best_hits(gene).itertuples():
                        other_seq_id = edge.seq1 if edge.seq1 != gene else edge.seq2
                        if other_seq_id.split("-")[0] == taxa:
                            other_seq_best_hits = self._get_best_hits(other_seq_id)
                            if gene in other_seq_best_hits.seq1.values or gene in other_seq_best_hits.seq2.values:
                                self.collapsed_genes.setdefault(gene, [])
                                self.collapsed_genes[gene].append(other_seq_id)
                                collapse_check_list.append(other_seq_id)
                                # Strip collapsed paralogs from the all-by-all graph
                                self.sim_scores = self.sim_scores[(self.sim_scores.seq1 != other_seq_id) &
                                                                  (self.sim_scores.seq2 != other_seq_id)]
                                if other_seq_id in collapsed_cluster:
                                    del collapsed_cluster[collapsed_cluster.index(other_seq_id)]
                                    del self.taxa[taxa][self.taxa[taxa].index(other_seq_id)]
                                    if other_seq_id in self.collapsed_genes:
                                        collapse_check_list += self.collapsed_genes[other_seq_id]
                                        del self.collapsed_genes[other_seq_id]
                                breakout = False
                                break
            self.seq_ids = collapsed_cluster

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
            raise ValueError("Parent of current cluster has not been named.")
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
            if gene not in self.seq_ids:
                raise NameError("Gene id '%s' not present in parental seq_ids" % gene)
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
                _valve = MyFuncs.SafetyValve(global_reps=10)
                while scores.score.std() == 0:
                    _valve.step("Failed to perturb:\n%s" % scores)
                    for _indx, _score in scores.score.iteritems():
                        scores.set_value(_indx, "score", random.gauss(_score, (_score * 0.0000001)))
                return scores

    def score(self):
        if self._score:
            return self._score
        # Confirm that a score for this cluster has not been caluclated before
        prev_scores = pd.read_csv(self.cluster_score_file.path, index_col=False)
        seq_ids = md5_hash("".join(sorted(self.seq_ids)))
        if seq_ids in prev_scores.cluster.values:
            self._score = float(prev_scores.score[prev_scores.cluster == seq_ids])
            return self._score

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

            # Get the similarity scores for within cliques and between cliques-remaining sequences, then generate kernel-density
            # functions for both. The overlap between the two functions is used to determine whether they should be separated
            self.cliques = []
            for clique in cliques:
                clique_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique)) & (self.sim_scores.seq2.isin(clique))]
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
                    clique = Cluster(clique, sim_scores=clique_scores, _parent=self, clique=True)
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

        self._score = self.raw_score(decliqued_cluster)
        for clique in self.cliques:
            if not clique:
                break
            clique_ids = md5_hash("".join(sorted(clique.seq_ids)))
            if clique_ids in prev_scores.cluster.values:
                self._score += float(prev_scores.score[prev_scores.cluster == clique_ids])
            else:
                clique_score = self.raw_score(clique.seq_ids)
                self._score += clique_score
                with self.lock:
                    self.cluster_score_file.write("\n%s,%s" % (clique_score, clique_ids))
                    sim_scores = self.sim_scores[(self.sim_scores.seq1.isin(clique.seq_ids)) &
                                                 (self.sim_scores.seq2.isin(clique.seq_ids))]
                    sim_scores.to_csv("%s/%s" % (self.similarity_graphs.path, clique_ids), index=False)
        with self.lock:
            self.cluster_score_file.write("\n%s,%s" % (self._score, seq_ids))
            self.sim_scores.to_csv("%s/%s" % (self.similarity_graphs.path, seq_ids), index=False)
        return self._score

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

        taxa = {}
        for gene in id_list:
            taxon = gene.split("-")[0]
            taxa.setdefault(taxon, [])
            taxa[taxon].append(gene)

        # An entire seq_ids should never consist of a single taxa because I've stripped out reciprocal best hit paralogs
        if len(taxa) == 1:
            raise ReferenceError("Only a single taxa found in seq_ids...")

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


def md5_hash(_input):
    _input = str(_input).encode()
    return md5(_input).hexdigest()


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


def _psi_pred(seq_obj, args):
    outdir = args[0]
    if os.path.isfile("%s/psi_pred/%s.ss2" % (outdir, seq_obj.id)):
        return
    temp_dir = MyFuncs.TempDir()
    pwd = os.getcwd()
    psipred_dir = os.path.abspath("%s/psipred" % os.path.dirname(__file__))
    os.chdir(temp_dir.path)
    with open("sequence.fa", "w") as _ofile:
        _ofile.write(seq_obj.format("fasta"))

    command = '''\
psiblast -db {0}/blastdb/pannexins -query sequence.fa -inclusion_ethresh 0.001 -out_pssm {1}/{2}.chk \
-num_iterations 3 -num_alignments 0 >& {1}/{2}.blast;
{0}/bin/chkparse {1}/{2}.chk > {1}/{2}.mtx;
{0}/bin/psipred {1}/{2}.mtx {0}/data/weights.dat {0}/data/weights.dat2 {0}/data/weights.dat3 > {1}/{2}.ss;
{0}/bin/psipass2 {0}/data/weights_p2.dat 1 1.0 1.0 {1}/{2}.ss2 {1}/{2}.ss > {1}/{2}.horiz;
'''.format(psipred_dir, temp_dir.path, seq_obj.id)

    Popen(command, shell=True).wait()
    os.chdir(pwd)
    shutil.move("%s/%s.ss2" % (temp_dir.path, seq_obj.id), "%s/psi_pred/" % outdir)
    return


def mcmcmc_mcl(args, params):
    inflation, gq = args
    external_tmp_dir, min_score, seqbuddy, parent_cluster = params
    mcl_tmp_dir = MyFuncs.TempDir()

    _output = Popen("mcl %s/input.csv --abc -te 2 -tf 'gq(%s)' -I %s -o %s/output.groups" %
                    (external_tmp_dir, gq, inflation, mcl_tmp_dir.path), shell=True, stderr=PIPE).communicate()

    _output = _output[1].decode()

    if re.search("\[mclvInflate\] warning", _output) and min_score:
        return min_score

    clusters = parse_mcl_clusters("%s/output.groups" % mcl_tmp_dir.path)
    score = 0
    prev_scores = pd.read_csv(parent_cluster.cluster_score_file.path, index_col=False)

    for indx, cluster in enumerate(clusters):
        cluster_ids = md5_hash("".join(sorted(cluster)))
        if cluster_ids in prev_scores.cluster.values:
            with parent_cluster.lock:
                sim_scores = pd.read_csv("%s/%s" % (parent_cluster.similarity_graphs.path, cluster_ids), index_col=False)
                scores = prev_scores.score[prev_scores.cluster == cluster_ids]
                if len(scores) > 1:
                    prev_scores = prev_scores.drop_duplicates()
                    prev_scores.to_csv(parent_cluster.cluster_score_file.path, index=False)
            score += float()
        else:
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % _id for _id in cluster]))
            alb_obj = generate_mas(sb_copy)
            sim_scores = create_all_by_all_scores(alb_obj, quiet=True)

        cluster = Cluster(cluster, sim_scores, _parent=parent_cluster)
        clusters[indx] = cluster
        score += cluster.score()

    with lock:
        with open("%s/max.txt" % external_tmp_dir, "r") as _ifile:
            current_max = float(_ifile.read())
        if score > current_max:
            write_mcl_clusters(clusters, "%s/best_group" % external_tmp_dir)
            with open("%s/max.txt" % external_tmp_dir, "w") as _ofile:
                _ofile.write(str(score))
    return score


def orthogroup_caller(master_cluster, cluster_list, seqbuddy, steps=1000, quiet=True):
    """
    Run MCMCMC on MCL to find the best orthogroups
    :param master_cluster: The group to be subdivided
    :type master_cluster: Cluster
    :param cluster_list: When a seq_ids is finalized after recursion, it is appended to this list
    :param seqbuddy: The sequences that are included in the master seq_ids
    :param steps: How many MCMCMC iterations to run TODO: calculate this on the fly
    :param quiet: Suppress StdErr
    :return: list of seq_ids objects
    """
    temp_dir = MyFuncs.TempDir()
    master_cluster.sim_scores.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    inflation_var = mcmcmc.Variable("I", 1.1, 20)
    gq_var = mcmcmc.Variable("gq", min(master_cluster.sim_scores.score), max(master_cluster.sim_scores.score))

    try:
        with open("%s/max.txt" % temp_dir.path, "w") as _ofile:
            _ofile.write("-1000000000")

        mcmcmc_factory = mcmcmc.MCMCMC([inflation_var, gq_var], mcmcmc_mcl, steps=steps, sample_rate=1,
                                       params=["%s" % temp_dir.path, False, seqbuddy, master_cluster], quiet=quiet,
                                       outfile="%s/mcmcmc_out.csv" % temp_dir.path)

    except RuntimeError:  # Happens when mcmcmc fails to find different initial chain parameters
        cluster_list.append(master_cluster)
        temp_dir.save("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name()))
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
        master_cluster.set_name()
        cluster_list.append(master_cluster)
        temp_dir.save("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name()))
        return cluster_list

    mcl_clusters = parse_mcl_clusters("%s/best_group" % temp_dir.path)
    for sub_cluster in mcl_clusters:
        cluster_ids = md5_hash("".join(sorted(sub_cluster)))
        sim_scores = pd.read_csv("%s/%s" % (master_cluster.similarity_graphs.path, cluster_ids), index_col=False)
        sub_cluster = Cluster(sub_cluster, sim_scores=sim_scores, _parent=master_cluster)
        if len(sub_cluster) in [1, 2]:
            sub_cluster.set_name()
            cluster_list.append(sub_cluster)
            continue

        seqbuddy_copy = Sb.make_copy(seqbuddy)
        seqbuddy_copy = Sb.pull_recs(seqbuddy_copy, ["^%s$" % rec_id for rec_id in sub_cluster.seq_ids])

        # Recursion...
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, steps=steps, quiet=quiet,)

    temp_dir.save("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name()))
    return cluster_list


def parse_mcl_clusters(path):
    with open(path, "r") as _ifile:
        clusters = _ifile.read()
    clusters = clusters.strip().split("\n")
    clusters = [cluster.strip().split("\t") for cluster in clusters]
    return clusters


def write_mcl_clusters(clusters, path):
    clusters_string = ""
    for cluster in clusters:
        clusters_string += "%s\n" % "\t".join(cluster.seq_ids)
    with open(path, "w") as _ofile:
        _ofile.write(clusters_string.strip())
    return


def merge_singles(clusters, scores):
    small_clusters = []
    small_group_names = []
    large_clusters = []
    large_group_names = []
    for cluster in clusters:
        if len(cluster.seq_ids) > 2:
            large_group_names.append(cluster.name())
            large_clusters.append(cluster)
        else:
            small_group_names.append(cluster.name())
            small_clusters.append(cluster)

    # Convert the large_clusters list to a dict using group name as key
    large_clusters = {x: large_clusters[j] for j, x in enumerate(large_group_names)}

    small_to_large_dict = {}
    for sclust in small_clusters:
        small_to_large_dict[sclust.name()] = {_ind: [] for _ind, value in large_clusters.items()}
        for sgene in sclust.seq_ids:
            for key, lclust in large_clusters.items():
                for lgene in lclust.seq_ids:
                    score = scores.loc[:][scores[0] == sgene]
                    score = score.loc[:][score[1] == lgene]

                    if score.empty:
                        score = scores.loc[:][scores[0] == lgene]
                        score = score.loc[:][score[1] == sgene]

                    if score.empty:
                        score = 0.
                    else:
                        score = float(score[2])

                    small_to_large_dict[sclust.name()][lclust.name()].append(score)

    small_clusters = {x: small_clusters[j] for j, x in enumerate(small_group_names)}
    for small_group_id, l_clusts in small_to_large_dict.items():
        # Convert data into list of numpy arrays that sm.stats can read, also get average scores for each seq_ids
        data = [np.array(x) for x in l_clusts.values()]
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
        result = sm.stats.multicomp.pairwise_tukeyhsd(df.observations, df.grouplabel)
        for line in str(result).split("\n")[4:-1]:
            line = re.sub("^ *", "", line.strip())
            line = re.sub(" +", "\t", line)
            line = line.split("\t")
            if max_ave in line:
                if line[5] == 'True':
                    continue
                else:
                    # Insufficient support to group the gene with max_ave group
                    break
        else:
            # The gene can be grouped with the max_ave group
            large_clusters[max_ave].seq_ids += small_clusters[small_group_id].seq_ids
            del small_clusters[small_group_id]

    clusters = [value for _ind, value in large_clusters.items()]
    clusters += [value for _ind, value in small_clusters.items()]
    return clusters
# NOTE: There used to be a support function. Check the GitHub history if there's desire to bring parts of it back


def score_sequences(_pair, args):
    # Calculate the best possible scores, and divide by the observed scores
    id1, id2 = _pair
    alb_obj, psi_pred_files, outfile = args
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
        with open(outfile, "a") as _ofile:
            _ofile.write("\n%s,%s,%s" % (id1, id2, final_score))
    return


def generate_mas(seqbuddy):
    if len(seqbuddy) == 1:
        alignment = Alb.AlignBuddy(str(seqbuddy))
    else:
        alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft", params="--globalpair --thread -1", quiet=True)
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

    # Don't want to modify the alignbuddy object in place
    alignment = Alb.make_copy(alignment)

    # Need to specify what columns the PsiPred files map to now that there are gaps.
    psi_pred_files = {}
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

    outfile = MyFuncs.TempFile()
    outfile.write("seq1,seq2,score")
    printer.clear()
    MyFuncs.run_multicore_function(all_by_all, score_sequences, [alignment, psi_pred_files, outfile.path], quiet=quiet)
    sim_scores = pd.read_csv(outfile.path, index_col=False)
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


def new_project():
    if os.path.exists(in_args.outdir):
        check = MyFuncs.ask("Output directory already exists, overwrite it [y]/n?") if not in_args.force else True
        if check:
            logging.info("Deleting all previous files from output directory.")
            shutil.rmtree(in_args.outdir)
            while os.path.exists(in_args.outdir):
                pass
        else:
            logging.warning("Program aborted by user to prevent overwriting of pre-existing output directory.")
            sys.exit()

    # Make all the directories needed for the run
    logging.info("mkdir %s" % in_args.outdir)
    os.makedirs(in_args.outdir)
    logging.info("mkdir %s/alignments" % in_args.outdir)
    os.makedirs("%s/alignments" % in_args.outdir)
    logging.info("mkdir %s/mcmcmc" % in_args.outdir)
    os.makedirs("%s/mcmcmc" % in_args.outdir)
    logging.info("mkdir %s/sim_scores" % in_args.outdir)
    os.makedirs("%s/sim_scores" % in_args.outdir)
    logging.info("mkdir %s/psi_pred\n" % in_args.outdir)
    os.makedirs("%s/psi_pred" % in_args.outdir)
    return


def resume():
    if not os.path.isdir(in_args.outdir):
        confirm = MyFuncs.ask("Specified input directory does not exist but the 'resume' flag was passed in. "
                              "Do you want to start a new whole new RD-MCL run [y]/n?: ")
        if confirm:
            new_project()
        else:
            logging.warning("Aborted... 'Resume' directory not found.")
            sys.exit()
    logging.info("RESUME: RD-MCL will attempt to load data.\n")
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
    parser.add_argument("-r", "--resume", help="Try to resume a job", action="store_true")

    in_args = parser.parse_args()

    logger_obj = Logger("rdmcl.log")
    logging.info("****************************** Recursive Dynamic Markov Clustering *******************************")
    logging.warning("RD-MCL version %s\n\n%s" % (VERSION, NOTICE))
    logging.info("**************************************************************************************************\n")
    logging.info("Working directory: %s" % os.getcwd())
    logging.info("Function call: %s\n" % " ".join(sys.argv))

    # Once passed this if/else, every step will try to read data from the output directory (i.e., attempt to 'resume').
    if in_args.resume:
        resume()
    else:
        new_project()

    # Move log file into output directory
    logger_obj.move_log("%s/rdmcl.log" % in_args.outdir)

    clique_check = True if not in_args.supress_clique_check else False
    recursion_check = True if not in_args.supress_recursion else False
    best = None
    best_clusters = None
    lock = Lock()
    printer = MyFuncs.DynamicPrint(quiet=in_args.quiet)
    timer = MyFuncs.Timer()

    sequences = Sb.SeqBuddy(in_args.sequences)
    BLOSUM62 = make_full_mat(SeqMat(MatrixInfo.blosum62))

    ambiguous_X = {"A": 0, "R": -1, "N": -1, "D": -1, "C": -2, "Q": -1, "E": -1, "G": -1, "H": -1, "I": -1, "L": -1,
                   "K": -1, "M": -1, "F": -1, "P": -2, "S": 0, "T": 0, "W": -2, "Y": -1, "V": -1}
    for aa in ambiguous_X:
        pair = sorted((aa, "X"))
        pair = tuple(pair)
        BLOSUM62[pair] = ambiguous_X[aa]

    # PSIPRED
    logging.warning("** PSI-Pred **")
    records_missing_ss_files = []
    records_with_ss_files = []
    for record in sequences.records:
        if os.path.isfile("%s/psi_pred/%s.ss2" % (in_args.outdir, record.id)):
            records_with_ss_files.append(record.id)
        else:
            records_missing_ss_files.append(record)
    if records_missing_ss_files and len(records_missing_ss_files) != len(sequences):
        logging.info("RESUME: PSI-Pred .ss2 files found for %s sequences:\n" % len(records_with_ss_files))

    if records_missing_ss_files:
        logging.warning("Executing PSI-Pred on %s sequences" % len(records_missing_ss_files))
        MyFuncs.run_multicore_function(records_missing_ss_files, _psi_pred, [in_args.outdir])
        logging.info("\tfinished in %s" % timer.split())
        logging.info("\tfiles saved to %s\n" % "%s/psi_pred/" % in_args.outdir)
    else:
        logging.warning("RESUME: All PSI-Pred .ss2 files found in %s/psi_pred/\n" % in_args.outdir)

    # Initial alignment
    logging.warning("** All-by-all graph **")
    gap_open = in_args.open_penalty
    gap_extend = in_args.extend_penalty
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (gap_open, gap_extend))

    if os.path.isfile("%s/alignments/group_0.aln" % in_args.outdir):
        logging.warning("RESUME: Initial multiple sequence alignment found")
        alignbuddy = Alb.AlignBuddy("%s/alignments/group_0.aln" % in_args.outdir)
    else:
        logging.warning("Generating initial multiple sequence alignment")
        alignbuddy = generate_mas(sequences)
        alignbuddy.write("%s/alignments/group_0.aln" % in_args.outdir)
        logging.info("\tfinished in %s" % timer.split())

    if os.path.isfile("%s/sim_scores/group_0.scores" % in_args.outdir):
        logging.warning("RESUME: Initial all-by-all similarity scores found")
        scores_data = pd.read_csv("%s/sim_scores/group_0.scores" % in_args.outdir, header=None,
                                  index_col=False, sep="\t")
        scores_data.columns = ['seq1', 'seq2', 'score']
    else:
        logging.warning("Generating initial all-by-all similarity graph")
        scores_data = create_all_by_all_scores(alignbuddy)
        scores_data.to_csv("%s/sim_scores/group_0.scores" % in_args.outdir, header=None, index=False, sep="\t")
        logging.info("\tfinished in %s" % timer.split())

    group_0 = pd.concat([scores_data.seq1, scores_data.seq2])
    group_0 = group_0.value_counts()
    group_0 = Cluster([i for i in group_0.index], scores_data)

    # Base cluster score
    logging.warning("** Scoring base cluster **")
    base_score = group_0.score()
    logging.warning("%s\n" % round(base_score, 4))

    #taxa_count = [x.split("-")[0] for x in master_cluster.seq_ids]
    #taxa_count = pd.Series(taxa_count)
    #taxa_count = taxa_count.value_counts()

    # Ortholog caller
    logging.warning("** Creating clusters **")
    final_clusters = []
    final_clusters = orthogroup_caller(group_0, final_clusters, seqbuddy=sequences,
                                       steps=in_args.mcmcmc_steps, quiet=False)

    logging.info("\tfinished in %s" % timer.split())

    # Fold singletons and doublets back into groups.
    if not in_args.supress_singlet_folding:
        logging.warning("** Folding singletons back into clusters **")
        final_clusters = merge_singles(final_clusters, scores_data)
        logging.warning("\tfinished in %s" % timer.split())

    # Format the clusters and output to stdout or file
    logging.warning("** Final formatting **")
    output = ""
    while len(final_clusters) > 0:
        _max = (0, 0)
        for ind, clust in enumerate(final_clusters):
            if len(clust.seq_ids) > _max[1]:
                _max = (ind, len(clust.seq_ids))

        ind, _max = _max[0], final_clusters[_max[0]]
        del final_clusters[ind]
        output += "group_%s\t%s\t" % (_max.name(), _max.score())
        for seq_id in _max.seq_ids:
            output += "%s\t" % seq_id
        output = "%s\n" % output.strip()

    logging.warning("\tfinished in %s\n" % timer.split())

    logging.warning("Total execution time: %s" % timer.total_elapsed())
    with open("%s/final_clusters.txt" % in_args.outdir, "w") as ofile:
        ofile.write(output)
        logging.warning("Final clusters written to: %s/final_clusters.txt" % in_args.outdir)
