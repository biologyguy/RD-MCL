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
from io import StringIO
from subprocess import Popen, PIPE, check_output, CalledProcessError
from multiprocessing import Lock, Pipe
from random import random
from math import ceil
from collections import OrderedDict

# 3rd party
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
from Bio.SubsMat import SeqMat, MatrixInfo

# My packages
import mcmcmc
import MyFuncs
import helpers
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb

# Globals
try:
    script_path = os.path.abspath(__file__).split("/")
    script_path = "/".join(script_path[:-1])
    git_commit = check_output(['git', '--git-dir=%s/../.git' % script_path, 'rev-parse', '--short', 'HEAD']).decode().strip()
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
CPUS = MyFuncs.usable_cpu_count()
TIMER = MyFuncs.Timer()


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
    def __init__(self, seq_ids, sim_scores, parent=None, clique=False):
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

        self.subgroup_counter = 0
        # self.clique_counter = 0
        # self.clique = clique
        # self.cliques = []
        self.cluster_score = None
        self.collapsed_genes = OrderedDict()  # If paralogs are reciprocal best hits, collapse them
        self._name = None
        for next_seq_id in seq_ids:
            taxa = next_seq_id.split(in_args.taxa_separator)[0]
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

            breakout = False
            while not breakout:
                breakout = True
                for seq1_id in seq_ids:
                    seq1_taxa = seq1_id.split(in_args.taxa_separator)[0]
                    paralog_best_hits = []
                    for hit in self._get_best_hits(seq1_id).itertuples():
                        seq2_id = hit.seq1 if hit.seq1 != seq1_id else hit.seq2
                        if seq2_id.split(in_args.taxa_separator)[0] != seq1_taxa:
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
        self.seq_id_hash = helpers.md5_hash("".join(sorted(self.seq_ids)))

    def _get_best_hits(self, gene):
        best_hits = self.sim_scores[(self.sim_scores.seq1 == gene) | (self.sim_scores.seq2 == gene)]
        if not best_hits.empty:
            best_hits = best_hits.loc[best_hits.score == max(best_hits.score)].values
            best_hits = pd.DataFrame(best_hits, columns=["seq1", "seq2", "score"])
        return best_hits

    def set_name(self):
        if self._name:
            pass
        elif not self.parent.name():
            raise ValueError("Parent of current cluster has not been named.\n%s" % self)
        # elif self.clique:
        #    self._name = "%s_c%s" % (self.parent.name(), self.parent.clique_counter)
        #    self.parent.clique_counter += 1
        else:
            self._name = "%s_%s" % (self.parent.name(), self.parent.subgroup_counter)
            self.parent.subgroup_counter += 1
            # for clique in self.cliques:
            #    clique.set_name()
        return

    def compare(self, query):
        matches = set(self.seq_ids).intersection(query.seq_ids)
        weighted_match = (len(matches) * 2.) / (len(self) + query.len)
        print("name: %s, matches: %s, weighted_match: %s" % (self.name(), len(matches), weighted_match))
        return weighted_match

    def name(self):
        if not self._name:
            raise AttributeError("Cluster has not been named.")
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
        self.cluster_score = self.raw_score(self.seq_ids)
        with LOCK:
            broker.query("""INSERT OR REPLACE INTO data_table (hash, graph, cluster_score, seq_ids, alignment)
  VALUES (  '{0}',
            '{1}',
            '{2}',
            (SELECT seq_ids FROM data_table WHERE hash = '{0}'),
            (SELECT alignment FROM data_table WHERE hash = '{0}')
          )""".format(self.seq_id_hash, self.sim_scores.to_csv(index=False), str(self.cluster_score)))
            # push(self.seq_id_hash, 'cluster_score', str(self.cluster_score))
            # push(self.seq_id_hash, 'graph', self.sim_scores.to_csv(index=False))
        return self.cluster_score

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

        if self.cliques:
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
            clique_ids = helpers.md5_hash("".join(sorted(clique.seq_ids)))
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in clique.seq_ids]))
            alb_obj = generate_msa(sb_copy)

            # if clique_ids in prev_scores:
            #    self.cluster_score += float(query(clique_ids, 'cluster_score'))
            # else:

            clique_score = self.raw_score(clique.seq_ids)
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

    @staticmethod
    def raw_score(id_list):
        if len(id_list) == 1:
            return 0

        taxa = OrderedDict()
        for gene in id_list:
            taxon = gene.split(in_args.taxa_separator)[0]
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
                unique_scores *= 2 * (1 + (len(group_0_cluster.taxa[taxon]) / len(group_0_cluster)))
            else:
                paralog_scores *= len(genes) * (len(genes) / len(group_0_cluster.taxa[taxon]))
        return unique_scores + paralog_scores

    def __len__(self):
        return len(self.seq_ids)

    def __str__(self):
        return str(self.seq_ids)


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
{0}/bin/seq2mtx sequence.fa > {1}/{2}.mtx;
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

    with LOCK:
        with open("%s/.progress" % in_args.outdir, "r") as ifile:
            _progress = json.load(ifile)
            _progress['mcl_runs'] += 1
        with open("%s/.progress" % in_args.outdir, "w") as ofile:
            json.dump(_progress, ofile)

    mcl_output = mcl_output[1].decode()
    if re.search("\[mclvInflate\] warning", mcl_output) and min_score:
        return min_score

    clusters = parse_mcl_clusters("%s/output.groups" % mcl_tmp_dir.path)
    score = 0

    for indx, cluster in enumerate(clusters):
        cluster_ids = helpers.md5_hash("".join(sorted(cluster)))
        cluster_data = broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(cluster_ids))
        if cluster_data and cluster_data[0][0]:
            with LOCK:
                sim_scores = pd.read_csv(StringIO(cluster_data[0][0]), index_col=False)
            score += float()
        else:
            sb_copy = Sb.make_copy(seqbuddy)
            sb_copy = Sb.pull_recs(sb_copy, "|".join(["^%s$" % rec_id for rec_id in cluster]))
            alb_obj = generate_msa(sb_copy)
            sim_scores = create_all_by_all_scores(alb_obj, quiet=True)

        cluster = Cluster(cluster, sim_scores, parent=parent_cluster)
        clusters[indx] = cluster
        score += cluster.score()

    with LOCK:
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
        cluster_list.append(master_cluster)
        if not os.path.isdir("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name())):
            temp_dir.save("%s/mcmcmc/%s" % (in_args.outdir, master_cluster.name()))
        alignment = generate_msa(seqbuddy)
        alignment.write("%s/alignments/%s.aln" % (in_args.outdir, master_cluster.name()))
        master_cluster.sim_scores.to_csv("%s/sim_scores/%s.scores" % (in_args.outdir, master_cluster.name()),
                                         header=None, index=False, sep="\t")

        with open("%s/.progress" % in_args.outdir, "r") as ifile:
            _progress = json.load(ifile)
            _progress["placed"] += len(master_cluster.seq_ids) if not master_cluster.subgroup_counter else 0
        with open("%s/.progress" % in_args.outdir, "w") as _ofile:
            json.dump(_progress, _ofile)
        return

    master_cluster.set_name()
    temp_dir = MyFuncs.TempDir()
    master_cluster.sim_scores.to_csv("%s/input.csv" % temp_dir.path, header=None, index=False, sep="\t")

    # If there are no paralogs in the cluster, then it is already at its highest score and MCL is unnecessary
    keep_going = False
    for taxon, genes in master_cluster.taxa.items():
        if len(genes) > 1:
            keep_going = True
            break
    if not keep_going:
        save_cluster()
        return cluster_list

    inflation_var = mcmcmc.Variable("I", 1.1, 20)
    gq_var = mcmcmc.Variable("gq", min(master_cluster.sim_scores.score), max(master_cluster.sim_scores.score))

    try:
        with open("%s/max.txt" % temp_dir.path, "w") as ofile:
            ofile.write("-1000000000")

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
        save_cluster()
        return cluster_list

    mcl_clusters = parse_mcl_clusters("%s/best_group" % temp_dir.path)
    recursion_clusters = []
    for sub_cluster in mcl_clusters:
        cluster_ids = helpers.md5_hash("".join(sorted(sub_cluster)))
        graph = broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(cluster_ids))[0][0]
        sim_scores = pd.read_csv(StringIO(graph), index_col=False)
        sub_cluster = Cluster(sub_cluster, sim_scores=sim_scores, parent=master_cluster)
        if sub_cluster.seq_id_hash == master_cluster.seq_id_hash:
            raise ArithmeticError("The sub_cluster and master_cluster are the same, but are returning different "
                                  "scores\nsub-cluster score: %s, master score: %s\n%s" % (sub_cluster.score(),
                                                                                           master_cluster.score(),
                                                                                           sub_cluster.seq_id_hash))
        sub_cluster.score()
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
        cluster_list = orthogroup_caller(sub_cluster, cluster_list, seqbuddy=seqbuddy_copy, steps=steps, quiet=quiet,)

    save_cluster()
    return cluster_list


def progress():
    with LOCK:
        with open("%s/.progress" % in_args.outdir, "r") as ifile:
            try:
                _progress = json.load(ifile)
                return "MCL runs processed: %s. Sequences placed: %s/%s. Run time: " % (_progress['mcl_runs'],
                                                                                        _progress['placed'],
                                                                                        _progress['total'])
            except json.decoder.JSONDecodeError:
                pass


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
            for taxon, seq_ids in small_clusters[sgroup].taxa.items():
                large_clusters[max_ave].taxa.setdefault(taxon, [])
                large_clusters[max_ave].taxa[taxon] += seq_ids
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


def score_sequences(seq_pairs, args):
    # Calculate the best possible scores, and divide by the observed scores
    alb_obj, psi_pred_files, output_dir = args
    file_name = helpers.md5_hash(str(seq_pairs))
    ofile = open("%s/%s" % (output_dir, file_name), "w")
    for seq_pair in seq_pairs:
        id1, id2 = seq_pair
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
        pair_score = (ss_score * 0.3) + (subs_mat_score * 0.7)
        ofile.write("\n%s,%s,%s" % (id1, id2, pair_score))
    ofile.close()
    return


def generate_msa(seqbuddy):
    seq_ids = sorted([rec.id for rec in seqbuddy.records])
    seq_id_hash = helpers.md5_hash("".join(seq_ids))
    alignment = broker.query("SELECT (alignment) FROM data_table WHERE hash='{0}'".format(seq_id_hash))
    if alignment and alignment[0][0]:
        alignment = Alb.AlignBuddy(alignment[0][0])
        pass
    else:
        if len(seqbuddy) == 1:
            alignment = Alb.AlignBuddy(str(seqbuddy))
        else:
            alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft", params="--globalpair --thread -2", quiet=True)

        broker.query("""INSERT OR REPLACE INTO data_table (hash, alignment, seq_ids, cluster_score, graph)
  VALUES (  '{0}',
            '{1}',
            '{2}',
            (SELECT cluster_score FROM data_table WHERE hash = '{0}'),
            (SELECT graph FROM data_table WHERE hash = '{0}')
          )""".format(seq_id_hash, str(alignment), str(", ".join(seq_ids))))

        # push(seq_id_hash, 'alignment', str(alignment))
        # push(seq_id_hash, 'seq_ids', str(", ".join(seq_ids)))
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
    seq_id_hash = helpers.md5_hash("".join(seq_ids))
    graph = broker.query("SELECT (graph) FROM data_table WHERE hash='{0}'".format(seq_id_hash))
    if graph and graph[0][0]:
        sim_scores = pd.read_csv(StringIO(graph[0][0]), index_col=False)
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

    all_by_all_outdir = MyFuncs.TempDir()
    if all_by_all:
        n = ceil(len(all_by_all) / CPUS)
        all_by_all = [all_by_all[i:i + n] for i in range(0, len(all_by_all), n)] if all_by_all else []
        MyFuncs.run_multicore_function(all_by_all, score_sequences, [alignment, psi_pred_files, all_by_all_outdir.path],
                                       quiet=quiet)
    sim_scores_file = MyFuncs.TempFile()
    sim_scores_file.write("seq1,seq2,score")
    aba_root, aba_dirs, aba_files = next(os.walk(all_by_all_outdir.path))
    for aba_file in aba_files:
        with open("%s/%s" % (aba_root, aba_file), "r") as ifile:
            sim_scores_file.write(ifile.read())
    sim_scores = pd.read_csv(sim_scores_file.get_handle("r"), index_col=False)
    return sim_scores


def check_sequences(seqbuddy):
    logging.warning("Checking that the format of all sequence ids matches 'taxa%sgene'" % in_args.taxa_separator)
    failures = []
    for rec in seqbuddy.records:
        rec_id = rec.id.split(in_args.taxa_separator)
        if len(rec_id) != 2:
            failures.append(rec.id)
    if failures:
        logging.error("Malformed sequence id(s): '%s'\nThe taxa separator character is currently set to '%s',\n"
                      " which can be changed with the '-ts' flag" % (", ".join(failures), in_args.taxa_separator))
        sys.exit()
    else:
        logging.warning("    %s sequences PASSED" % len(seqbuddy))
    return

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(prog="orthogroup_caller", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("sequences", help="Location of sequence file", action="store")
    parser.add_argument("outdir", action="store", default="%s/rd-mcd" % os.getcwd(),
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
                        type=float, default=-5)
    parser.add_argument("-ep", "--extend_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                        type=float, default=0)
    parser.add_argument("-ts", "--taxa_separator", action="store", default="-",
                        help="Specify a string that separates taxa ids from gene names")
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
    logging.info("Working directory: %s" % os.getcwd())
    if not shutil.which("mcl"):
        logging.error("The 'mcl' program is not detected on your system (see http://micans.org/mcl/).")
        sys.exit()
    mcl = Popen("mcl --version", stdout=PIPE, shell=True).communicate()[0].decode()
    logging.info("MCL version: %s" % re.search("mcl (.*)", mcl).group(1))

    if not shutil.which("mafft"):
        logging.error("The 'MAFFT' program is not detected "
                      "on your system (see http://mafft.cbrc.jp/alignment/software/).")
        sys.exit()
    mafft = Popen("mafft --version", stderr=PIPE, shell=True).communicate()[1].decode()
    logging.info("MAFFT version: %s" % mafft.strip())

    if not os.path.isdir(in_args.outdir):
        logging.info("Creating output directory: %s" % in_args.outdir)
        os.makedirs(in_args.outdir)

    logging.info("")  # Just want to insert an extra line break for asthetics
    logging.warning("Launching SQLite Daemon")

    broker = helpers.SQLiteBroker(db_file="{0}/sqlite_db.sqlite".format(in_args.outdir))
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()
    sequences = Sb.SeqBuddy(in_args.sequences)
    check_sequences(sequences)
    seq_ids_hash = helpers.md5_hash("".join(sorted([rec.id for rec in sequences.records])))

    # Check if job has been run already
    if broker.query("SELECT (hash) FROM data_table WHERE hash='{0}'".format(seq_ids_hash)) == seq_ids_hash:
        logging.warning("RESUME: This output directory was previous used for an identical RD-MCL run.\n"
                        "        All cached resources will be reused.")

    # Make sure all the necessary directories are present and emptied of old run files
    for outdir in ["%s%s" % (in_args.outdir, x) for x in ["", "/alignments", "/mcmcmc", "/sim_scores", "/psi_pred"]]:
        if not os.path.isdir(outdir):
            logging.info("mkdir %s" % outdir)
            os.makedirs(outdir)
        # Delete old 'group' files/directories
        root, dirs, files = next(os.walk(outdir))
        for _file in files:
            if "group" in _file:
                os.remove("%s/%s" % (root, _file))
        for _dir in dirs:
            if "group" in _dir:
                shutil.rmtree("%s/%s" % (root, _dir))

    # Move log file into output directory
    logger_obj.move_log("%s/rdmcl.log" % in_args.outdir)

    clique_check = True if not in_args.supress_clique_check else False
    recursion_check = True if not in_args.supress_recursion else False

    BLOSUM62 = helpers.make_full_mat(SeqMat(MatrixInfo.blosum62))

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
        logging.info("\t-- finished in %s --" % TIMER.split())
        logging.info("\tfiles saved to %s" % "%s/psi_pred/" % in_args.outdir)
    else:
        logging.warning("RESUME: All PSI-Pred .ss2 files found in %s/psi_pred/" % in_args.outdir)

    # Initial alignment
    logging.warning("\n** All-by-all graph **")
    gap_open = in_args.open_penalty
    gap_extend = in_args.extend_penalty
    logging.info("gap open penalty: %s\ngap extend penalty: %s" % (gap_open, gap_extend))

    align_data = broker.query("SELECT (alignment) FROM data_table WHERE hash='{0}'".format(seq_ids_hash))
    if align_data and align_data[0][0]:
        logging.warning("RESUME: Initial multiple sequence alignment found")
        alignbuddy = Alb.AlignBuddy(align_data[0][0])
    else:
        logging.warning("Generating initial multiple sequence alignment with MAFFT")
        alignbuddy = generate_msa(sequences)
        alignbuddy.write("%s/alignments/group_0.aln" % in_args.outdir)

        logging.info("\t-- finished in %s --" % TIMER.split())
    if os.path.isfile("%s/sim_scores/complete_all_by_all.scores" % in_args.outdir):
        logging.warning("RESUME: Initial all-by-all similarity graph found")
        scores_data = pd.read_csv("%s/sim_scores/complete_all_by_all.scores" % in_args.outdir,
                                  index_col=False, sep="\t")
        scores_data.columns = ["seq1", "seq2", "score"]
        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data)

    else:
        logging.warning("Generating initial all-by-all similarity graph")
        logging.info(" written to: %s/sim_scores/complete_all_by_all.scores" % in_args.outdir)
        scores_data = create_all_by_all_scores(alignbuddy)
        scores_data.to_csv("%s/sim_scores/complete_all_by_all.scores" % in_args.outdir,
                           header=None, index=False, sep="\t")
        broker.query("UPDATE data_table SET graph='{0}' "
                     "WHERE hash='{1}'".format(scores_data.to_csv(header=None, index=False), seq_ids_hash))

        group_0_cluster = Cluster([rec.id for rec in sequences.records], scores_data)
        logging.info("\t-- finished in %s --" % TIMER.split())

    # Base cluster score
    base_score = group_0_cluster.score()
    logging.warning("Base cluster score: %s" % round(base_score, 4))
    # Reset the score of group_0 to account for possible collapse of paralogs
    group_0_cluster.score()
    if group_0_cluster.collapsed_genes:
        logging.warning("Reciprocal best-hit cliques of paralogs have been identified in the input sequences.")
        logging.info(" A representative sequence has been selected from each clique, and the remaining")
        logging.info(" sequences will be placed in the final clusters at the end of the run.")
        with open("%s/paralog_cliques.json" % in_args.outdir, "w") as outfile:
            json.dump(group_0_cluster.collapsed_genes, outfile)
            logging.warning(" Cliques written to: %s/paralog_cliques.json" % in_args.outdir)

    # taxa_count = [x.split(in_args.taxa_separator)[0] for x in group_0_cluster.seq_ids]
    # taxa_count = pd.Series(taxa_count)
    # taxa_count = taxa_count.value_counts()

    # Ortholog caller
    logging.warning("\n** Recursive MCL **")
    final_clusters = []
    with open("%s/.progress" % in_args.outdir, "w") as progress_file:
        progress_dict = {"mcl_runs": 0, "placed": 0, "total": len(group_0_cluster)}
        json.dump(progress_dict, progress_file)
    run_time = MyFuncs.RunTime(prefix=progress, _sleep=0.3, final_clear=True)
    run_time.start()
    final_clusters = orthogroup_caller(group_0_cluster, final_clusters, seqbuddy=sequences,
                                       steps=in_args.mcmcmc_steps, quiet=True)
    run_time.end()

    with open("%s/.progress" % in_args.outdir, "r") as progress_file:
        progress_dict = json.load(progress_file)
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
    # Fold singletons and doublets back into groups. This can't be 'resumed', because it changes the clusters
    if not in_args.supress_singlet_folding:
        logging.warning("\n** Folding orphan sequences into clusters **")
        final_clusters = place_orphans(final_clusters, group_0_cluster.sim_scores)
        logging.warning("\t-- finished in %s --" % TIMER.split())

    # Format the clusters and output to file
    logging.warning("\n** Final formatting **")
    if group_0_cluster.collapsed_genes:
        logging.warning("Placing collapsed paralogs into their respective clusters")
        for clust in final_clusters:
            for gene_id, paralogs in group_0_cluster.collapsed_genes.items():
                if gene_id in clust.seq_ids:
                    clust.seq_ids += paralogs

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
    with open("%s/final_clusters.txt" % in_args.outdir, "w") as outfile:
        outfile.write(output)
        logging.warning("Final score: %s" % round(final_score, 4))
        if final_score < base_score:
            logging.warning("    This is lower than the initial base score after"
                            " the inclusion of collapsed paralogs")
        logging.warning("Clusters written to: %s/final_clusters.txt" % in_args.outdir)

    broker.close()
