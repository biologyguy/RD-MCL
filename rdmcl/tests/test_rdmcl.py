#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from .. import rdmcl
from .. import helpers
import os
import sqlite3
import pandas as pd
from collections import OrderedDict
from buddysuite import buddy_resources as br
from copy import copy, deepcopy
import shutil
import argparse


# #########  Mock classes and functions  ########## #
class MockLogging(object):
    @staticmethod
    def warning(_input):
        print(_input)

    @staticmethod
    def error(_input):
        print(_input)


# #########  Cluster class and functions  ########## #
def test_cluster_instantiate_group_0(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    assert [taxa for taxa in cluster.taxa] == ['BOL', 'Bab', 'Bch', 'Bfo', 'Bfr', 'Cfu', 'Dgl', 'Edu', 'Hca', 'Hru',
                                               'Hvu', 'Lcr', 'Lla', 'Mle', 'Oma', 'Pba', 'Tin', 'Vpa']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "ae11f765826e65eeefb52c3d04776d9b"
    assert cluster.taxa_sep == "-"
    assert cluster.parent is None
    assert cluster.subgroup_counter == 0
    assert cluster.cluster_score is None
    assert cluster.collapsed_genes == OrderedDict([('Lcr-PanxαA', ['Lcr-PanxαL']),
                                                   ('Mle-Panxα10A', ['Mle-Panxα9']),
                                                   ('Hvu-PanxβL', ['Hvu-PanxβO', 'Hvu-PanxβA', 'Hvu-PanxβC',
                                                                   'Hvu-PanxβB', 'Hvu-PanxβJ', 'Hvu-PanxβM',
                                                                   'Hvu-PanxβI', 'Hvu-PanxβF', 'Hvu-PanxβH',
                                                                   'Hvu-PanxβD', 'Hvu-PanxβG', 'Hvu-PanxβE',
                                                                   'Hvu-PanxβK'])])
    assert cluster._name == "group_0"
    assert cluster.seq_ids == ['BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαE', 'BOL-PanxαF',
                               'BOL-PanxαG', 'BOL-PanxαH', 'Bab-PanxαA', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαD',
                               'Bab-PanxαE', 'Bch-PanxαA', 'Bch-PanxαB', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαE',
                               'Bfo-PanxαA', 'Bfo-PanxαB', 'Bfo-PanxαC', 'Bfo-PanxαD', 'Bfo-PanxαE', 'Bfo-PanxαF',
                               'Bfo-PanxαG', 'Bfo-PanxαH', 'Bfo-PanxαI', 'Bfo-PanxαJ', 'Bfr-PanxαA', 'Bfr-PanxαB',
                               'Bfr-PanxαC', 'Bfr-PanxαD', 'Cfu-PanxαA', 'Cfu-PanxαB', 'Cfu-PanxαC', 'Cfu-PanxαD',
                               'Cfu-PanxαE', 'Cfu-PanxαF', 'Dgl-PanxαA', 'Dgl-PanxαB', 'Dgl-PanxαC', 'Dgl-PanxαD',
                               'Dgl-PanxαE', 'Dgl-PanxαF', 'Dgl-PanxαG', 'Dgl-PanxαH', 'Dgl-PanxαI', 'Edu-PanxαA',
                               'Edu-PanxαB', 'Edu-PanxαC', 'Edu-PanxαD', 'Edu-PanxαE', 'Edu-PanxαF', 'Edu-PanxαG',
                               'Edu-PanxαH', 'Edu-PanxαI', 'Hca-PanxαA', 'Hca-PanxαB', 'Hca-PanxαC', 'Hca-PanxαD',
                               'Hca-PanxαE', 'Hca-PanxαF', 'Hca-PanxαG', 'Hca-PanxαH', 'Hru-PanxαA', 'Hru-PanxαB',
                               'Hru-PanxαC', 'Hru-PanxαD', 'Hru-PanxαE', 'Hvu-PanxβL', 'Lcr-PanxαA', 'Lcr-PanxαB',
                               'Lcr-PanxαC', 'Lcr-PanxαD', 'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαG', 'Lcr-PanxαH',
                               'Lcr-PanxαI', 'Lcr-PanxαJ', 'Lcr-PanxαK', 'Lla-PanxαA', 'Lla-PanxαB', 'Lla-PanxαC',
                               'Mle-Panxα1', 'Mle-Panxα10A', 'Mle-Panxα11', 'Mle-Panxα12', 'Mle-Panxα2', 'Mle-Panxα3',
                               'Mle-Panxα4', 'Mle-Panxα5', 'Mle-Panxα6', 'Mle-Panxα7A', 'Mle-Panxα8', 'Oma-PanxαA',
                               'Oma-PanxαB', 'Oma-PanxαC', 'Oma-PanxαD', 'Pba-PanxαA', 'Pba-PanxαB', 'Pba-PanxαC',
                               'Pba-PanxαD', 'Pba-PanxαE', 'Pba-PanxαF', 'Pba-PanxαG', 'Tin-PanxαA', 'Tin-PanxαB',
                               'Tin-PanxαC', 'Tin-PanxαD', 'Tin-PanxαE', 'Tin-PanxαF', 'Vpa-PanxαA', 'Vpa-PanxαB',
                               'Vpa-PanxαC', 'Vpa-PanxαD', 'Vpa-PanxαE', 'Vpa-PanxαF', 'Vpa-PanxαG']
    assert cluster.seq_id_hash == "0eb00400b766199af2040597b963f4da"

    with pytest.raises(ValueError) as err:
        sim_scores = hf.get_data("cteno_sim_scores")
        sim_scores = sim_scores.ix[1:, :]
        rdmcl.Cluster(cluster.seq_ids, sim_scores)
    assert "The number of incoming sequence ids (119) does not match the expected graph size of 7021" in str(err)


def test_cluster_instantiate_child(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    sim_scores = pd.read_csv("%sCteno_pannexins_subgroup_sim.scores" % hf.resource_path, index_col=False, header=None)
    sim_scores.columns = ["seq1", "seq2", "score"]
    cluster = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert [taxa for taxa in cluster.taxa] == ['BOL', 'Bab', 'Bch', 'Bfo', 'Dgl', 'Edu', 'Hca',
                                               'Hru', 'Lcr', 'Mle', 'Oma', 'Tin', 'Vpa']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "eb1ef296e9f4e74cd6deea490a447326"
    assert cluster.taxa_sep == "-"
    assert cluster.parent._name == "group_0"
    assert cluster.subgroup_counter == 0
    assert cluster.cluster_score is None
    assert cluster.collapsed_genes == OrderedDict([('Mle-Panxα10A', ['Mle-Panxα9'])])
    assert cluster._name is None
    assert cluster.seq_ids == child_ids
    assert cluster.seq_id_hash == hf.string2hash(", ".join(sorted(child_ids)))


def test_cluster_get_name(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    assert parent.name() == "group_0"
    sub_cluster_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC']
    child = rdmcl.Cluster(sub_cluster_ids, hf.get_sim_scores(sub_cluster_ids), parent=parent)
    with pytest.raises(AttributeError) as err:
        child.name()
    assert "Cluster has not been named." in str(err)


def test_cluster_set_name(hf):
    group_0 = rdmcl.Cluster(*hf.base_cluster_args())
    assert group_0._name == "group_0"
    group_0.set_name()
    assert group_0._name == "group_0"

    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    sim_scores = pd.read_csv("%sCteno_pannexins_subgroup_sim.scores" % hf.resource_path, index_col=False, header=None)
    sim_scores.columns = ["seq1", "seq2", "score"]
    group_0_0 = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=group_0)
    assert group_0_0._name is None

    grandchild_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB']
    group_0_0_0 = rdmcl.Cluster(grandchild_ids, hf.get_sim_scores(grandchild_ids), parent=group_0_0)
    with pytest.raises(ValueError) as err:
        group_0_0_0.set_name()
    assert "Parent of current cluster has not been named." in str(err)

    group_0_0.set_name()
    assert group_0_0._name == "group_0_0"


def test_cluster_compare(hf, capsys):
    subject_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE']
    subject = rdmcl.Cluster(subject_ids, hf.get_sim_scores(subject_ids))

    query_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    query = rdmcl.Cluster(query_ids, hf.get_sim_scores(query_ids))
    assert subject.compare(query) == 0.2
    out, err = capsys.readouterr()
    assert out == "name: group_0, matches: 1, weighted_match: 0.2\n"
    assert query.compare(subject) == 0.2


def test_cluster_get_best_hits(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    best_hit = cluster.get_best_hits("Bab-PanxαA")
    assert best_hit.iloc[0].seq2 == "Lcr-PanxαG"


def test_cluster_recursive_best_hits(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    global_best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
    best_hits = cluster.recursive_best_hits('Bab-PanxαB', global_best_hits, ['Bab-PanxαB'])
    assert best_hits.to_csv() == """\
,psi,raw_score,score,seq1,seq2,subsmat
0,0.9841450231425388,0.9710485816574068,0.9776353490763704,Bab-PanxαB,Vpa-PanxαB,0.9748454887622982
1,0.9832934130111736,0.98033753188254,0.9867858512031742,Lcr-PanxαH,Vpa-PanxαB,0.9882826104283174
2,0.9832934130111736,0.98033753188254,0.9867858512031742,Lcr-PanxαH,Vpa-PanxαB,0.9882826104283174
""", print(best_hits.to_csv())


def test_cluster_perturb(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    global_best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
    best_hits = cluster.recursive_best_hits('Lcr-PanxαH', global_best_hits, ['Lcr-PanxαH'])
    assert best_hits.iloc[0].score == 0.9867858512031742

    for indx, score in best_hits["score"].iteritems():
        best_hits.set_value(indx, "score", 0.9867858512031742)

    best_hits = cluster.perturb(best_hits)
    assert best_hits.iloc[0].score != 0.9867858512031742
    assert round(best_hits.iloc[0].score, 5) == 0.98679


def test_cluster_get_base_cluster(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)

    # No paralogs
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    # First time calling Cluster.score() calculates the score
    assert child.score() == 46.138159371492705

    # The second call just retrieves the attribute from the cluster saved during first call
    assert child.score() == 46.138159371492705

    # With paralogs
    child_ids = ['BOL-PanxαA', 'BOL-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child.score() == 41.36858164983165

    # Single sequence
    child_ids = ['BOL-PanxαA']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert round(child.score(), 12) == 1.847222222222

    # Include an orphan sequence
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC']
    sim_scores = hf.get_sim_scores(child_ids)
    child = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert round(child.score(), 3) == 8.575
    child.seq_ids.append("Foo-Bar3")
    assert round(child.score(force=True), 12) == 8.510526315789

    # Edge case where child is full size of parent
    child = rdmcl.Cluster(parent.seq_ids, parent.sim_scores, parent=parent)
    assert round(child.score(), 12) == 264.57969700078


def test_cluster_pull_scores_subgraph(hf):
    clust = rdmcl.Cluster(*hf.base_cluster_args())
    output = str(clust.pull_scores_subgraph(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC']))
    assert output == """\
            seq1        seq2   subsmat       psi  raw_score     score
132   Hru-PanxαA  Lcr-PanxαH  0.974896  0.952851   0.965774  0.968282
165   Hru-PanxαA  Tin-PanxαC  0.968037  0.937704   0.958389  0.958937
6851  Lcr-PanxαH  Tin-PanxαC  0.984287  0.968946   0.975098  0.979684""", print(output)


def test_rbhc_less_than_6(hf):
    # Cluster of less than 6 seqs returns self
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_0 ####### #
['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC']
\tTERMINATED: Group too small."""


def test_rbhc_all_seqs_in_clique(hf):
    # If clique brings in entire cluster, return self
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    seq_ids = ['Mle-Panxα10A', 'Mle-Panxα9', 'Vpa-PanxαB', 'BOL-PanxαA', 'Dgl-PanxαE', 'Lcr-PanxαH']
    sim_scores = hf.get_sim_scores(seq_ids)

    # Modify the scores to create a complete clique
    indx = sim_scores[(sim_scores.seq2 == 'Mle-Panxα10A') & (sim_scores.seq1 == 'Mle-Panxα9')]
    sim_scores.set_value(indx.index[0], "raw_score", 0.7)

    indx = sim_scores[(sim_scores.seq1 == 'Mle-Panxα10A') & (sim_scores.seq2 == 'Vpa-PanxαB')]
    sim_scores.set_value(indx.index[0], "raw_score", 0.8)

    indx = sim_scores[(sim_scores.seq2 == 'Vpa-PanxαB') & (sim_scores.seq1 == 'BOL-PanxαA')]
    sim_scores.set_value(indx.index[0], "raw_score", 0.99)

    indx = sim_scores[(sim_scores.seq2 == 'Mle-Panxα9') & (sim_scores.seq1 == 'BOL-PanxαA')]
    sim_scores.set_value(indx.index[0], "raw_score", 0.99)

    indx = sim_scores[(sim_scores.seq1 == 'Mle-Panxα9') & (sim_scores.seq2 == 'Vpa-PanxαB')]
    sim_scores.set_value(indx.index[0], "raw_score", 0.97)

    log_file = br.TempFile()
    cluster = rdmcl.Cluster(seq_ids, sim_scores, parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_0 ####### #
['BOL-PanxαA', 'Dgl-PanxαE', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα9', 'Vpa-PanxαB']
\tTERMINATED: Entire cluster pulled into clique on Mle.""", print(log_file.read().strip())


def test_rbhc_no_paraogs(hf):
    # If no paralogs in cluster, return self
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE', 'Vpa-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_0 ####### #
['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαE']
\tTERMINATED: No paralogs present."""


def test_rbhc_all_disqualified(hf):
    # When all cliques are disqualified due to size and/or overlap
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
               'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα9', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_0 ####### #
['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA',\
 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα9', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
\tChecking taxa for overlapping cliques:
\t\t# #### Mle #### #
\t\tDisqualified cliques:
\t\t\t['Mle-Panxα9', 'Vpa-PanxαB'] is too small
\t\t\t['Mle-Panxα10A', 'Mle-Panxα9', 'Vpa-PanxαB'] overlaps with ['Mle-Panxα9', 'Vpa-PanxαB']

\t\t!! ALL CLIQUES DISQUALIFIED !!

\tTERMINATED: No significant cliques identified."""


def test_rbhc_creates_orphan(hf):
    # When spinning off an otherwise valid clique would orphan one or more sequences, return self
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI',
               'Bch-PanxαB',  'Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert "TERMINATED: Spinning off cliques would orphan Bch-PanxαB and Lcr-PanxαI" in log_file.read()


def test_rbhc_fail_integration(hf):
    # If the outer scope sim scores are not sufficiently separated from inner scope
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI', 'Lcr-PanxαL', 'Mle-Panxα11']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    cliques = cluster.create_rbh_cliques(log_file)
    assert cliques[0] == cluster
    assert """Test for KDE separation:
\t\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Lcr-PanxαL', 'Mle-Panxα11', 'Mle-Panxα4']
\t\t\tOuter KDE: {'shape': (1, 5), 'covariance': 0.0320703176923, 'inv_cov': 31.1814809443,\
 '_norm_factor': 2.2444584476}
\t\t\tClique KDE: {'shape': (1, 10), 'covariance': 0.0196599709781, 'inv_cov': 50.8647749845,\
  '_norm_factor': 3.51464423219}""" in log_file.read()
    assert "FAIL" in log_file.read()


def test_rbhc_multi_clique_pass(hf):
    # All sequences spread across valid cliques
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI', 'Lcr-PanxαL', 'Mle-Panxα5']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    cliques = cluster.create_rbh_cliques(log_file)
    assert cliques[0].seq_ids == ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4']
    assert cliques[1].seq_ids == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαA']

    assert """Test for KDE separation:
\t\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4']
\t\t\tOuter KDE: {'shape': (1, 9), 'covariance': 0.0247415314604, 'inv_cov': 40.4178699123,\
 '_norm_factor': 3.54850754302}
\t\t\tClique KDE: {'shape': (1, 3), 'covariance': 9.93711384681e-05, 'inv_cov': 10063.2841227,\
  '_norm_factor': 0.0749620270178}
""" in log_file.read(), print(log_file.read())

    assert """Cliques identified and spun off:
\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4']
\t\t['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαA']""" in log_file.read()


def test_rbhc_combine_cliques_pass(hf):
    # Check whether identified cliques can be combined
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI', 'Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    cliques = cluster.create_rbh_cliques(log_file)
    assert cliques[0].seq_ids == ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4', 'Vpa-PanxαA']
    assert cliques[1].seq_ids == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    assert """Checking if new cliques can combine
\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4'] - ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']  NO
\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4'] - ['BOL-PanxαF', 'Mle-Panxα4', 'Vpa-PanxαA']  YES""" in log_file.read()


def test_rbhc_remaining_seqs_pass(hf):
    # Cliques identified, but there are some remaining sequences that need to be bundled together
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()
    seq_ids = ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF', 'Edu-PanxαG', 'Pba-PanxαE', 'Lcr-PanxαA']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    new_clusters = cluster.create_rbh_cliques(log_file)
    assert new_clusters[0].seq_ids == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    assert new_clusters[1].seq_ids == ['Edu-PanxαG', 'Lcr-PanxαA', 'Pba-PanxαE']
    assert """\t\t\tPASS

\tCliques identified and spun off:
\t\t['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
\t\t['Edu-PanxαG', 'Lcr-PanxαA', 'Pba-PanxαE']""" in log_file.read()


def test_cluster_len(hf):
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids))
    assert len(cluster) == 5


def test_cluster_str(hf):
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids))
    assert str(cluster) == "['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC']"


def test_cluster2database(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()
    rdmcl.cluster2database(cluster, broker, ">Seq1\nMPQQCS-SS\n>Seq2\nMPQICMAAS")
    broker.close()

    connect = sqlite3.connect("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    cursor = connect.cursor()
    cursor.execute("SELECT * FROM data_table")
    response = cursor.fetchall()
    assert len(response) == 1
    assert response[0][0] == "0eb00400b766199af2040597b963f4da"                 # hash
    assert 'BOL-PanxαA, BOL-PanxαB, BOL-PanxαC, BOL-PanxαD' in response[0][1]   # seq_ids
    assert response[0][2] == '>Seq1\nMPQQCS-SS\n>Seq2\nMPQICMAAS'               # alignment
    assert 'Hca-PanxαG,Lla-PanxαC,0.291912106162889' in response[0][3]          # graph
    assert response[0][4] == '264.57969700077973'                               # score
    connect.close()


# #########  PSI-PRED  ########## #
def test_run_psi_pred(hf):
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Bab-PanxαA")
    ss2_file = rdmcl.run_psi_pred(seqbuddy.records[0])
    assert hf.string2hash(ss2_file) == "b50f39dc22e4d16be325efdd14f7900d"


def test_mc_psi_pred(hf):
    tmpdir = br.TempDir()
    rdmcl.mc_psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open("%s/BOL-PanxαB.ss2" % tmpdir.path, "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.mc_psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "af9666d37426caa2bbf6b9075ce8df96", print(output)


def test_read_ss2_file(hf):
    ss2 = rdmcl.read_ss2_file("%spsi_pred%sMle-Panxα10A.ss2" % (hf.resource_path, hf.sep))
    assert type(ss2) == pd.DataFrame
    assert list(ss2.columns) == ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
    assert len(ss2.index) == 429
    assert str(ss2.loc[1]) == """\
indx              2
aa                R
ss                C
coil_prob     0.906
helix_prob    0.037
sheet_prob    0.028
Name: 1, dtype: object"""


def test_compare_psi_pred(hf):
    ss2_1 = rdmcl.read_ss2_file("%spsi_pred%sMle-Panxα10A.ss2" % (hf.resource_path, hf.sep))
    ss2_2 = rdmcl.read_ss2_file("%spsi_pred%sMle-Panxα8.ss2" % (hf.resource_path, hf.sep))
    assert rdmcl.compare_psi_pred(ss2_1, ss2_2) == 0.691672882672883


# #########  Orthogroup caller  ########## #
def test_orthogroup_caller(hf):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    seq_ids = ['BOL-PanxαA', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαF', 'BOL-PanxαG', 'BOL-PanxαH', 'Lcr-PanxαA',
               'Lcr-PanxαB', 'Lcr-PanxαD', 'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαH', 'Lcr-PanxαI', 'Lcr-PanxαJ',
               'Lcr-PanxαK', 'Lcr-PanxαL', 'Mle-Panxα1', 'Mle-Panxα2', 'Mle-Panxα4', 'Mle-Panxα5', 'Mle-Panxα6',
               'Mle-Panxα7A', 'Mle-Panxα8', 'Mle-Panxα9', 'Mle-Panxα10A', 'Mle-Panxα12', 'Vpa-PanxαA', 'Vpa-PanxαB',
               'Vpa-PanxαC', 'Vpa-PanxαF', 'Vpa-PanxαG']

    seqbuddy = rdmcl.Sb.SeqBuddy(hf.resource_path + "Cteno_pannexins.fa")
    rdmcl.Sb.pull_recs(seqbuddy, "^%s$" % "$|^".join(seq_ids))

    cluster = rdmcl.Cluster([rec.id for rec in seqbuddy.records],
                            hf.get_db_graph("ec5f4340bd57eb9db6819802598457c7", broker), collapse=True)

    cluster_list = []
    outdir = br.TempDir()
    outdir.subdir("progress")
    outdir.subdir("alignments")
    outdir.subdir("mcmcmc")
    outdir.subdir("psi_pred")
    outdir.subdir("sim_scores")

    progress = rdmcl.Progress(os.path.join(outdir.path, "progress"), cluster)
    psi_pred_ss2_dfs = OrderedDict()
    for rec in seqbuddy.records:
        psi_pred_ss2_dfs[rec.id] = rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, os.sep, rec.id))
    steps = 10
    r_seed = 2
    orthogroups = rdmcl.orthogroup_caller(cluster, cluster_list, seqbuddy, broker, progress, outdir.path,
                                          psi_pred_ss2_dfs, steps=steps, r_seed=r_seed)

    expected = [['BOL-PanxαC', 'Lcr-PanxαE', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα6'],
                ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4', 'Vpa-PanxαA'],
                ['Lcr-PanxαA', 'Mle-Panxα5', 'Vpa-PanxαF'],
                ['BOL-PanxαC', 'Lcr-PanxαH', 'Mle-Panxα10A'],
                ['BOL-PanxαD', 'Lcr-PanxαD', 'Mle-Panxα2'],
                ['Lcr-PanxαK', 'Mle-Panxα7A'],
                ['Lcr-PanxαB', 'Mle-Panxα1'],
                ['BOL-PanxαH', 'Mle-Panxα8'],
                ['BOL-PanxαG', 'Lcr-PanxαF'],
                ['Vpa-PanxαC'],
                ['Mle-Panxα6'],
                ['Lcr-PanxαE']]

    orthogroup_seqs = [clust.seq_ids for clust in orthogroups]
    assert cluster.seq_ids in orthogroup_seqs
    for clust in expected:
        assert clust in orthogroup_seqs


# #########  Miscellaneous  ########## #
def test_progress(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args(), collapse=True)
    tmpdir = br.TempDir()
    progress = rdmcl.Progress(tmpdir.path, cluster)
    assert progress.outdir == tmpdir.path
    assert os.path.isfile("{0}{1}.progress".format(tmpdir.path, hf.sep))
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # The dictionary is not static, so just sort the string:
        # {"placed": 0, "mcl_runs": 0, "total": 119}
        assert "".join(sorted(ifile.read())) == '     """""",,00119:::_aaccdelllmnoprsttu{}'

    progress.update("mcl_runs", 2)
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # {"placed": 0, "mcl_runs": 2, "total": 119}
        assert "".join(sorted(ifile.read())) == '     """""",,01129:::_aaccdelllmnoprsttu{}'

    json = progress.read()
    assert json["mcl_runs"] == 2
    assert json["placed"] == 0
    assert json["total"] == 119

    assert str(progress) == "MCL runs processed: 2. Sequences placed: 0/119. Run time: "


def test_check_sequences(hf, monkeypatch, capsys):
    monkeypatch.setattr(rdmcl, "logging", MockLogging)
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.check_sequences(seqbuddy, "-")
    out, err = capsys.readouterr()
    assert "Checking that the format of all sequence ids matches 'taxa-gene'" in out
    assert "    134 sequences PASSED" in out

    rdmcl.Sb.rename(seqbuddy, "Mle-Panxα1", "Mle:Panxα1")
    assert not rdmcl.check_sequences(seqbuddy, "-")

    out, err = capsys.readouterr()
    assert "Checking that the format of all sequence ids matches 'taxa-gene'" in out
    print(out)
    assert "Malformed sequence id(s): 'Mle:Panxα1, Mle:Panxα10A, Mle:Panxα11, Mle:Panxα12'\n" \
           "The taxa separator character is currently set to '-',\n" \
           " which can be changed with the '-ts' flag" in out

"""
def test_generate_msa(hf):
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()

    broker.query("INSERT INTO data_table (hash) VALUES ('991d38af45b2b71022eb6348679db953')")
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    all_mle_alignment = rdmcl.generate_msa(seqbuddy, broker)

    broker.query("INSERT INTO data_table (hash) VALUES ('c2bfd6a538e8876c5770a8f07b6b220e')")
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "9")
    mle9_alignment = rdmcl.generate_msa(seqbuddy, broker)

    connect = sqlite3.connect("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    cursor = connect.cursor()
    cursor.execute("SELECT * FROM data_table")
    response = cursor.fetchall()
    for indx in range(5):  # There's race condition here on the database
        if len(response) != 2 or not (response[1][2]):
            helpers.sleep(0.5)
            cursor.execute("SELECT * FROM data_table")
            response = cursor.fetchall()
        else:
            break

    assert len(response) == 2
    assert hf.string2hash(response[0][2]) == "22ba0f62bb616d1106f0a43ac73d343e"
    assert hf.string2hash(str(all_mle_alignment)) == "22ba0f62bb616d1106f0a43ac73d343e"

    assert hf.string2hash(response[1][2]) == '919d26d7db868d01fa285090eb98299e'
    assert hf.string2hash(str(mle9_alignment)) == "919d26d7db868d01fa285090eb98299e"

    alignment = rdmcl.generate_msa(seqbuddy, broker)
    assert hf.string2hash(str(alignment)) == "919d26d7db868d01fa285090eb98299e"

    cursor.execute("SELECT * FROM data_table")
    response = cursor.fetchall()
    assert len(response) == 2
    connect.close()
    broker.stop_broker()
"""


def test_mc_score_sequence(hf):
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    alignbuddy = rdmcl.Alb.generate_msa(seqbuddy, "mafft", params="--globalpair --thread 2", quiet=True)
    psi_pred_files = [(rec.id, rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, hf.sep, rec.id)))
                      for rec in alignbuddy.records()]
    psi_pred_files = OrderedDict(psi_pred_files)
    tmp_dir = br.TempDir()

    seq_pairs = [("Mle-Panxα2", "Mle-Panxα12"), ("Mle-Panxα1", "Mle-Panxα3")]
    args = [alignbuddy, psi_pred_files, tmp_dir.path, -5, 0]
    rdmcl.mc_score_sequences(seq_pairs, args)
    assert os.path.isfile("%s%s50687872cdaaed1fee7e809b5a032377" % (tmp_dir.path, hf.sep))
    with open("%s%s50687872cdaaed1fee7e809b5a032377" % (tmp_dir.path, hf.sep), "r") as ifile:
        content = ifile.read()
        assert content == """
Mle-Panxα2,Mle-Panxα12,0.3760737547928561,0.6767437070938215
Mle-Panxα1,Mle-Panxα3,0.2237341940045855,0.6656420029895364""", print(content)


def test_compare_pairwise_alignment():
    seqbuddy = rdmcl.Sb.SeqBuddy(">seq1\nMPQMSASWI\n>Seq2\nMPPQISASI")
    alignbuddy = rdmcl.Alb.generate_msa(seqbuddy, "mafft", params="--globalpair --op 0 --thread 2", quiet=True)
    assert str(alignbuddy) == """\
>seq1
MP-QMSASWI
>Seq2
MPPQISAS-I
"""
    subs_mat_score = rdmcl.compare_pairwise_alignment(alignbuddy, -5, 0)
    assert subs_mat_score == 0.4658627087198515


def test_create_all_by_all_scores(hf):
    tmpdir = br.TempDir()
    sql_broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    sql_broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                           "graph TEXT", "cluster_score TEXT"])
    sql_broker.start_broker()

    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    psi_pred_files = [(rec.id, rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, hf.sep, rec.id)))
                      for rec in seqbuddy.records]
    psi_pred_files = OrderedDict(psi_pred_files)

    rdmcl.ALIGNMETHOD = "mafft"
    rdmcl.ALIGNPARAMS = "--globalpair"
    sim_scores, alignbuddy = rdmcl.create_all_by_all_scores(seqbuddy, psi_pred_files, sql_broker)
    assert len(sim_scores.index) == 66  # This is for 12 starting sequences --> (a * (a - 1)) / 2
    compare = sim_scores.loc[:][(sim_scores['seq1'] == "Mle-Panxα2") & (sim_scores['seq2'] == "Mle-Panxα12")]
    assert """\
          seq1         seq2  subsmat      psi  raw_score    score
42  Mle-Panxα2  Mle-Panxα12  0.31444  0.25294    0.50618  0.29599""" in str(compare)
    assert len(alignbuddy.records()) == 12

    rdmcl.Sb.pull_recs(seqbuddy, "Mle-Panxα2")
    sim_scores, alignbuddy = rdmcl.create_all_by_all_scores(seqbuddy, psi_pred_files, sql_broker)
    assert len(sim_scores.index) == 0
    sql_broker.close()


# #########  MCL stuff  ########## #
def test_mcmcmc_mcl(hf):
    ext_tmp_dir = br.TempDir()
    ext_tmp_dir.subfile("max.txt")
    cluster_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Hca-PanxαB',
                   'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Tin-PanxαC',
                   'Vpa-PanxαB', 'Oma-PanxαC', 'Edu-PanxαA', 'Bch-PanxαC']
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    seqbuddy = rdmcl.Sb.pull_recs(seqbuddy, "^%s$" % "$|^".join(cluster_ids))
    taxa_sep = "-"
    sql_broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    sql_broker.start_broker()
    
    cluster = rdmcl.Cluster(cluster_ids, hf.get_db_graph("3c15516819aa19b069b0e8858444f876", sql_broker), collapse=True)
    psi_pred_ss2_dfs = OrderedDict()
    for rec in seqbuddy.records:
        psi_pred_ss2_dfs[rec.id] = rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, os.sep, rec.id))
    os.makedirs(os.path.join(ext_tmp_dir.path, "progress"))
    progress = rdmcl.Progress(os.path.join(ext_tmp_dir.path, "progress"), cluster)

    args = (6.372011782427792, 0.901221218627, 1)  # inflation, gq, r_seed
    params = [ext_tmp_dir.path, seqbuddy, cluster, taxa_sep, sql_broker, psi_pred_ss2_dfs, progress, 3]

    assert rdmcl.mcmcmc_mcl(args, params) == 18.307692307692307
    with open(os.path.join(ext_tmp_dir.path, "max.txt"), "r") as ifile:
        output = ifile.read()
        assert output == "76d0cf67bbdd773b2676a2bab2e0c1c8,cb23cf3b4d355140e525a1158af5102d," \
                         "8521872b6c07205f3198bb70699f3d93,a06adee8cc3631773890bb5842bf8df9," \
                         "3dee08bfba869d1128cedd1680dfefe6,09cac4f034df8a2805171e1e61cc8666", print(output)

    args = (3.1232, 0.73432, 1)  # inflation, gq, r_seed
    assert rdmcl.mcmcmc_mcl(args, params) == 19.538461538461537
    with open(os.path.join(ext_tmp_dir.path, "max.txt"), "r") as ifile:
        output = ifile.read()
        assert output == "6b39ebc4f5fe7dfef786d8ee3e1594ed,cb23cf3b4d355140e525a1158af5102d," \
                         "8521872b6c07205f3198bb70699f3d93,a06adee8cc3631773890bb5842bf8df9," \
                         "09cac4f034df8a2805171e1e61cc8666\n" \
                         "76d0cf67bbdd773b2676a2bab2e0c1c8,cb23cf3b4d355140e525a1158af5102d," \
                         "8521872b6c07205f3198bb70699f3d93,a06adee8cc3631773890bb5842bf8df9," \
                         "3dee08bfba869d1128cedd1680dfefe6,09cac4f034df8a2805171e1e61cc8666", print(output)

    args = (10.1232, 0.43432, 1)  # inflation, gq, r_seed
    assert rdmcl.mcmcmc_mcl(args, params) == 24.153846153846153
    with open(os.path.join(ext_tmp_dir.path, "max.txt"), "r") as ifile:
        output = ifile.read()
        assert output == "", print(output)

    with open(os.path.join(ext_tmp_dir.path, "best_group"), "r") as ifile:
        output = ifile.read()
        assert output == "BOL-PanxαA	Bab-PanxαB	Bch-PanxαC	Bfo-PanxαB	Dgl-PanxαE	Hca-PanxαB	Hru-PanxαA	" \
                         "Lcr-PanxαH	Mle-Panxα10A	Oma-PanxαC	Tin-PanxαC	Vpa-PanxαB\n" \
                         "Edu-PanxαA", print(output)
    sql_broker.close()


def test_parse_mcl_clusters(hf):
    clusters = rdmcl.parse_mcl_clusters("%sCteno_pannexins_mcl_clusters.clus" % hf.resource_path)
    assert clusters[6] == ["BOL-PanxαH", "Dgl-PanxαH", "Edu-PanxαC", "Hca-PanxαF", "Mle-Panxα8", "Pba-PanxαC"]


def test_write_mcl_clusters(hf):
    clusters = rdmcl.parse_mcl_clusters("%sCteno_pannexins_mcl_clusters.clus" % hf.resource_path)
    clusters = [rdmcl.Cluster(cluster, hf.get_sim_scores(cluster)) for cluster in clusters]
    tmp_file = br.TempFile()
    rdmcl.write_mcl_clusters(clusters, tmp_file.path)
    assert hf.string2hash(tmp_file.read()) == "4024a3aa15d605fd8924d8fc6cb9361e"


# #########  Orphan placement  ########## #
def test_instantiate_orphan(hf):
    # Get starting graph from pre-computed database
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    psi_pred_ss2_dfs = hf.get_data("ss2_dfs")

    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    graph = hf.get_db_graph("6935966a6b9967c3006785488d968230", broker)

    parent_ids = [rec.id for rec in parent_sb.records]
    parent_cluster = rdmcl.Cluster(parent_ids, graph, collapse=False)

    cluster1 = ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA']
    # alignbuddy = rdmcl.generate_msa(rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb),
    #                                                   "^%s$" % "$|^".join(cluster1)), broker)
    graph = hf.get_db_graph("14f1cd0e985ed87b4e31bc07453481d2", broker)
    # graph = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_ss2_dfs)
    cluster1 = rdmcl.Cluster(cluster1, graph, parent=parent_cluster)
    # rdmcl.cluster2database(cluster1, broker, alignbuddy)
    cluster1.set_name()

    cluster2 = ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD', 'Pba-PanxαB', 'Tin-PanxαF']
    # alignbuddy = rdmcl.generate_msa(rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb),
    #                                                   "^%s$" % "$|^".join(cluster2)), broker)
    graph = hf.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker)
    # graph = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_ss2_dfs)
    cluster2 = rdmcl.Cluster(cluster2, graph, parent=parent_cluster)
    # rdmcl.cluster2database(cluster2, broker, alignbuddy)
    cluster2.set_name()

    cluster3 = ['Vpa-PanxαD']  # This should be placed in cluster 2
    # alignbuddy = rdmcl.generate_msa(rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb),
    #                                                   "^%s$" % "$|^".join(cluster3)), broker)
    graph = hf.get_db_graph("ded0bc087974589c24945bb36197d36f", broker)
    # graph = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_ss2_dfs)
    cluster3 = rdmcl.Cluster(cluster3, graph, parent=parent_cluster)
    # rdmcl.cluster2database(cluster3, broker, alignbuddy)
    cluster3.set_name()

    cluster4 = ['Hca-PanxαA', 'Lcr-PanxαG']  # This should be placed in cluster 1
    # alignbuddy = rdmcl.generate_msa(rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb),
    #                                                   "^%s$" % "$|^".join(cluster4)), broker)
    graph = hf.get_db_graph("035b3933770942c7e32e27c06e619825", broker)
    # graph = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_ss2_dfs)
    cluster4 = rdmcl.Cluster(cluster4, graph, parent=parent_cluster)
    # rdmcl.cluster2database(cluster4, broker, alignbuddy)
    cluster4.set_name()

    cluster5 = ['Hvu-PanxβA']  # This should not be placed in a cluster
    # alignbuddy = rdmcl.generate_msa(rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb),
    #                                                   "^%s$" % "$|^".join(cluster5)), broker)
    graph = hf.get_db_graph("c62b6378d326c2479296c98f0f620d0f", broker)
    # graph = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_ss2_dfs)
    cluster5 = rdmcl.Cluster(cluster5, graph, parent=parent_cluster)
    # rdmcl.cluster2database(cluster5, broker, alignbuddy)
    cluster5.set_name()
    clusters = [cluster1, cluster2, cluster3, cluster4, cluster5]

    broker.close()

    # Have Orphans create the subgraph from scratch
    tmpdir = br.TempDir()
    broker = helpers.SQLiteBroker("%s%sdb.sqlite" % (tmpdir.path, hf.sep))
    broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                       "graph TEXT", "cluster_score TEXT"])
    broker.start_broker()
    orphans = rdmcl.Orphans(parent_sb, clusters, broker, psi_pred_ss2_dfs)

    assert orphans.seqbuddy == parent_sb
    assert orphans.clusters == clusters
    assert orphans.sql_broker == broker
    assert orphans.psi_pred_ss2_dfs == psi_pred_ss2_dfs
    assert orphans.num_orphans == 0
    small_clusters = [clust.seq_ids for clust_name, clust in orphans.small_clusters.items()]
    assert small_clusters == [cluster3.seq_ids, cluster4.seq_ids, cluster5.seq_ids]
    large_clusters = [clust.seq_ids for clust_name, clust in orphans.large_clusters.items()]
    assert large_clusters == [cluster1.seq_ids, cluster2.seq_ids]
    assert len([seq_id for sub_clust in small_clusters + large_clusters for seq_id in sub_clust]) == 14
    assert type(orphans.tmp_file) == rdmcl.br.TempFile
    all_sim_scores_str = str(orphans.lrg_cluster_sim_scores)
    assert round(orphans.lrg_cluster_sim_scores.iloc[0], 12) == 0.961941661285
    assert round(orphans.lrg_cluster_sim_scores.iloc[-1], 12) == 0.940557899149
    assert len(orphans.lrg_cluster_sim_scores) == 20  # This is for the two large clusters --> Σ(a*(a-1))/2
    broker.close()

    # Run a second time, reading the graph from the open database this time (happens automatically)
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    orphans = rdmcl.Orphans(seqbuddy, clusters, broker, psi_pred_ss2_dfs)
    assert all_sim_scores_str == str(orphans.lrg_cluster_sim_scores)
    assert len(orphans.lrg_cluster_sim_scores) == 20
    broker.close()

# ToDo: Broken!
def test_place_orphans(hf):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()

    graph = hf.get_db_graph("6935966a6b9967c3006785488d968230", broker)
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)

    psi_pred_ss2_dfs = OrderedDict()
    for rec in parent_sb.records:
        psi_pred_ss2_dfs[rec.id] = rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, os.sep, rec.id))

    parent_ids = [rec.id for rec in parent_sb.records]
    parent_cluster = rdmcl.Cluster(parent_ids, graph, collapse=False)

    cluster1 = ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA']
    cluster1 = rdmcl.Cluster(cluster1, hf.get_db_graph("14f1cd0e985ed87b4e31bc07453481d2", broker),
                             parent=parent_cluster)
    cluster1.set_name()

    cluster2 = ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD', 'Pba-PanxαB', 'Tin-PanxαF']
    cluster2 = rdmcl.Cluster(cluster2, hf.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker),
                             parent=parent_cluster)
    cluster2.set_name()

    cluster3 = ['Vpa-PanxαD']  # This should be placed in cluster 2
    cluster3 = rdmcl.Cluster(cluster3, hf.get_db_graph("ded0bc087974589c24945bb36197d36f", broker),
                             parent=parent_cluster)
    cluster3.set_name()

    cluster4 = ['Hca-PanxαA', 'Lcr-PanxαG']  # This should be placed in cluster 1
    cluster4 = rdmcl.Cluster(cluster4, hf.get_db_graph("035b3933770942c7e32e27c06e619825", broker),
                             parent=parent_cluster)
    cluster4.set_name()

    cluster5 = ['Hvu-PanxβA']  # This should not be placed in a cluster
    cluster5 = rdmcl.Cluster(cluster5, hf.get_db_graph("c62b6378d326c2479296c98f0f620d0f", broker),
                             parent=parent_cluster)
    cluster5.set_name()
    orig_clusters = [cluster1, cluster2, cluster3, cluster4, cluster5]
    clusters = [deepcopy(clust) for clust in orig_clusters]

    orphans = rdmcl.Orphans(parent_sb, clusters, broker, psi_pred_ss2_dfs)

    # ##### mc_check_orphans() ##### #
    tmp_file = br.TempFile()
    orphans.mc_check_orphans(cluster3, [tmp_file.path])
    orphans.mc_check_orphans(cluster4, [tmp_file.path])
    orphans.mc_check_orphans(cluster5, [tmp_file.path])
    assert tmp_file.read() == """\
group_0_2\tgroup_0_1\t0.0374077704349
group_0_3\tgroup_0_0\t0.0192173134297
""", print(tmp_file.read())

    assert hf.string2hash(orphans.tmp_file.read()) == "e10539a2921508a1c8c8707ee14b4409", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### _check_orphan() ##### #
    # Cluster is placed
    orphs = orphans._check_orphan(cluster3)

    assert orphs[0] == 'group_0_1'
    assert round(orphs[1], 12) == 0.037407770435
    assert not orphans._check_orphan(cluster5)  # Insufficient support for largest cluster
    assert hf.string2hash(orphans.tmp_file.read()) == "a9646c7748aacbf3f0553fdae26a23e9", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### place_orphans() ##### #
    orphans.place_orphans(multi_core=False)
    assert orphans.clusters[0].seq_ids == ["Hvu-PanxβA"]
    assert orphans.clusters[1].seq_ids == ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE',
                                           'Bfr-PanxαA', 'Hca-PanxαA', 'Lcr-PanxαG']
    assert orphans.clusters[2].seq_ids == ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD',
                                           'Pba-PanxαB', 'Tin-PanxαF', 'Vpa-PanxαD']
    assert hf.string2hash(orphans.tmp_file.read()) == "3535657d136a774c9c2d0161f5d2dcd6", print(orphans.tmp_file.read())

    # Multicore doesn't seem to work in py.test, but can at least call it like it does
    clusters = [deepcopy(clust) for clust in orig_clusters]
    orphans = rdmcl.Orphans(parent_sb, clusters, broker, psi_pred_ss2_dfs)
    orphans.place_orphans(multi_core=True)
    assert len(orphans.clusters) == 3

    # If no small clusters, nothing much happens
    orphans = rdmcl.Orphans(parent_sb, clusters[:2], broker, psi_pred_ss2_dfs)
    assert not orphans.place_orphans(multi_core=False)
    assert len(orphans.clusters) == 2

    # ToDo: Test for failure on Tukey HSD test


# #########  User Interface  ########## #
parser = argparse.ArgumentParser(prog="orthogroup_caller", description="",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("sequences", help="Location of sequence file", action="store", nargs="*")
parser.add_argument("outdir", action="store", default=os.path.join(os.getcwd(), "rd-mcd"), nargs="*",
                    help="Where should results be written?")
parser.add_argument("-sql", "--sqlite_db", action="store", help="Specify a SQLite database location.")
parser.add_argument("-psi", "--psipred_dir", action="store",
                    help="If PSI-Pred files are pre-calculated, tell us where.")
parser.add_argument("-mcs", "--mcmc_steps", default=1000, type=int,
                    help="Specify how deeply to sample MCL parameters")
parser.add_argument("-sr", "--suppress_recursion", action="store_true",
                    help="Stop after a single round of MCL. For testing.")
parser.add_argument("-scc", "--suppress_clique_check", action="store_true",
                    help="Do not check for or break up cliques. For testing.")
parser.add_argument("-ssf", "--suppress_singlet_folding", action="store_true",
                    help="Do not check for or merge singlets. For testing.")
parser.add_argument("-sit", "--suppress_iteration", action="store_true",
                    help="Only check for cliques and orphans once")
parser.add_argument("-spc", "--suppress_paralog_collapse", action="store_true",
                    help="Do not merge best hit paralogs")
parser.add_argument("-nt", "--no_msa_trim", action="store_true",
                    help="Don't apply the gappyout algorithm to MSAs before scoring")
parser.add_argument("-op", "--open_penalty", help="Penalty for opening a gap in pairwise alignment scoring",
                    type=float, default=rdmcl.GAP_OPEN)
parser.add_argument("-ep", "--ext_penalty", help="Penalty for extending a gap in pairwise alignment scoring",
                    type=float, default=rdmcl.GAP_EXTEND)
parser.add_argument("-ts", "--taxa_sep", action="store", default="-",
                    help="Specify a string that separates taxa ids from gene names")
parser.add_argument("-cpu", "--max_cpus", type=int, action="store", default=rdmcl.CPUS,
                    help="Specify the maximum number of cores RD-MCL can use.")
parser.add_argument("-drb", "--dr_base", help="Set the base for diminishing returns function", type=float)
parser.add_argument("-cnv", "--converge", type=float,
                    help="Set minimum Gelman-Rubin PSRF value for convergence")
parser.add_argument("-rs", "--r_seed", help="Specify a random seed for repeating a specific run", type=int)
parser.add_argument("-ch", "--chains", default=rdmcl.MCMC_CHAINS, type=int,
                    help="Specify how many MCMCMC chains to run (default=3)")
parser.add_argument("-wlk", "--walkers", default=3, type=int,
                    help="Specify how many Metropolis-Hastings walkers are in each chain (default=2)")
parser.add_argument("-lwt", "--lock_wait_time", type=int, default=1200, metavar="",
                    help="Specify num seconds a process should wait on the SQLite database before crashing"
                         " out (default=1200)")
parser.add_argument("-wdb", "--workdb", action="store", default="",
                    help="Specify the location of a sqlite database monitored by workers")
parser.add_argument("-algn_m", "--align_method", action="store", default="clustalo",
                          help="Specify which alignment algorithm to use (supply full path if not in $PATH)")
parser.add_argument("-algn_p", "--align_params", action="store", default="",
                          help="Supply alignment specific parameters")
parser.add_argument("-f", "--force", action="store_true",
                    help="Overwrite previous run")
parser.add_argument("-q", "--quiet", action="store_true",
                    help="Suppress all output during run (only final output is returned)")

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


def test_argparse_init(monkeypatch, hf):
    out_dir = br.TempDir()
    argv = ['rdmcl.py', os.path.join(hf.resource_path, "BOL_Lcr_Mle_Vpa.fa"), out_dir.path]
    monkeypatch.setattr(rdmcl.sys, "argv", argv)
    temp_in_args = rdmcl.argparse_init()
    assert temp_in_args.mcmc_steps == 0
    assert temp_in_args.open_penalty == -5
    assert temp_in_args.ext_penalty == 0


def test_full_run(hf, capsys):
    # I can't break these up into separate test functions because of collisions with logger
    out_dir = br.TempDir()

    '''
    # First try a RESUME run
    test_in_args = deepcopy(in_args)
    test_in_args.sequences = os.path.join(hf.resource_path, "BOL_Lcr_Mle_Vpa.fa")
    test_in_args.outdir = out_dir.path
    test_in_args.sqlite_db = os.path.join(hf.resource_path, "db.sqlite")
    test_in_args.psipred_dir = os.path.join(hf.resource_path, "psi_pred")
    test_in_args.mcmc_steps = 10
    test_in_args.r_seed = 1
    rdmcl.full_run(test_in_args)

    for expected_dir in ["alignments", "mcmcmc", "psi_pred", "sim_scores"]:
        assert os.path.isdir(os.path.join(out_dir.path, expected_dir))

    for expected_file in ["cliques.log", "final_clusters.txt", "orphans.log", "paralog_cliques", "rdmcl.log"]:
        assert os.path.isfile(os.path.join(out_dir.path, expected_file))

    with open(os.path.join(out_dir.path, "final_clusters.txt"), "r") as ifile:
        content = ifile.read()
        assert content == """\
group_0_0_0\t22.5308\tBOL-PanxαA\tBOL-PanxαC\tLcr-PanxαE\tLcr-PanxαH\tLcr-PanxαJ\tMle-Panxα10A\tMle-Panxα12\t\
Mle-Panxα6\tMle-Panxα9\tVpa-PanxαB\tVpa-PanxαC\tVpa-PanxαG
group_0_1\t15.6667\tBOL-PanxαF\tLcr-PanxαI\tMle-Panxα4\tVpa-PanxαA
group_0_2\t11.0069\tLcr-PanxαA\tLcr-PanxαL\tMle-Panxα5\tVpa-PanxαF
group_0_3\t7.875\tBOL-PanxαD\tLcr-PanxαD\tMle-Panxα2
group_0_4\t3.75\tLcr-PanxαB\tMle-Panxα1
group_0_5\t4.875\tBOL-PanxαH\tMle-Panxα8
group_0_6\t4.875\tBOL-PanxαG\tLcr-PanxαF
group_0_0_1\t3.75\tLcr-PanxαK\tMle-Panxα7A
""", print(content)

    out, err = capsys.readouterr()
    assert "RESUME: All PSI-Pred .ss2 files found" in err
    assert "RESUME: Initial multiple sequence alignment found" in err
    assert "RESUME: Initial all-by-all similarity graph found" in err
    assert "Iterative placement of orphans and paralog RBHC removal" in err
    '''
    # Now a full run from scratch (on smaller set), with non-existant psi-pred dir
    open(os.path.join(out_dir.path, "rdmcl.log"), "w").close()
    shutil.move(os.path.join(out_dir.path, "rdmcl.log"), "rdmcl.log")
    test_in_args = deepcopy(in_args)
    out_dir = br.TempDir()
    # Make some files and a directory that will need to be removed
    out_dir.subdir("mcmcmc")
    out_dir.subdir(os.path.join("mcmcmc", "group_0"))
    out_dir.subdir("psi_pred")
    out_dir.subfile("group_0.txt")
    out_dir.subfile("orphans.log")
    out_dir.subfile("cliques.log")
    subset_ids = ["BOL-PanxαA", "Lcr-PanxαH", "Mle-Panxα10A", "Mle-Panxα9", "Vpa-PanxαB", "BOL-PanxαF",
                  "Lcr-PanxαI", "Mle-Panxα4", "Vpa-PanxαA", "BOL-PanxαC", "Mle-Panxα12", "Vpa-PanxαG"]
    for seq in subset_ids[:-1]:
        shutil.copyfile(os.path.join(hf.resource_path, "psi_pred", seq + ".ss2"),
                        os.path.join(out_dir.path, "psi_pred", seq + ".ss2"))
    seqbuddy = rdmcl.Sb.SeqBuddy(os.path.join(hf.resource_path, "BOL_Lcr_Mle_Vpa.fa"))
    rdmcl.Sb.pull_recs(seqbuddy, "^%s$" % "$|^".join(subset_ids))
    seqbuddy.write(os.path.join(out_dir.path, "seqbuddy"))
    test_in_args.sequences = os.path.join(out_dir.path, "seqbuddy")
    test_in_args.outdir = out_dir.path
    test_in_args.psipred_dir = "foo-bared_psi_pred"  # This doesn't exist
    test_in_args.mcmc_steps = 10
    test_in_args.r_seed = 1
    rdmcl.full_run(test_in_args)

    for expected_dir in ["alignments", "mcmcmc", "psi_pred", "sim_scores"]:
        assert os.path.isdir(os.path.join(out_dir.path, expected_dir))

    for expected_file in ["cliques.log", "final_clusters.txt", "orphans.log",
                          "paralog_cliques", "rdmcl.log", "sqlite_db.sqlite"]:
        assert os.path.isfile(os.path.join(out_dir.path, expected_file))

    with open(os.path.join(out_dir.path, "final_clusters.txt"), "r") as ifile:
        content = ifile.read()
        assert content == """\
group_0_0\t12.2708\tBOL-PanxαA\tLcr-PanxαH\tMle-Panxα10A\tMle-Panxα9\tVpa-PanxαB
group_0_1\t12.0\tBOL-PanxαF\tLcr-PanxαI\tMle-Panxα4\tVpa-PanxαA
group_0_2\t7.0\tBOL-PanxαC\tMle-Panxα12\tVpa-PanxαG
""", print(content)

    out, err = capsys.readouterr()
    assert "Generating initial all-by-all similarity graph (66 comparisons)" in err
