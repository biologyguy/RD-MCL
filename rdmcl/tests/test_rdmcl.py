#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import sqlite3
import pandas as pd
import shutil
import argparse
import numpy as np
from .. import rdmcl
from .. import helpers
from math import ceil
from collections import OrderedDict
from buddysuite import buddy_resources as br
from copy import deepcopy

pd.set_option('expand_frame_repr', False)


# #########  Mock classes and functions  ########## #
class MockLogging(object):
    @staticmethod
    def warning(_input):
        print(_input)

    @staticmethod
    def error(_input):
        print(_input)

    @staticmethod
    def info(_input):
        print(_input)


class MockPopen(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def wait(self):
        print(self.args, self.kwargs)
        return


class MockHeartBeat(object):
    def __init__(self, hbdb_path, pulse_rate, thread_type='master', dummy=False):
        self.hbdb_path = hbdb_path
        self.pulse_rate = pulse_rate
        self.id = None
        self.running_process = None
        self.check_file = br.TempFile()
        self.thread_type = thread_type
        self.dummy = dummy

    def start(self):
        return

    def end(self):
        return


def mock_keyboardinterupt(*args, **kwargs):
    raise KeyboardInterrupt(args, kwargs)


# #########  Cluster class and functions  ########## #
def test_cluster_instantiate_group_0(hf, monkeypatch):
    monkeypatch.setattr(rdmcl.Cluster, "collapse", lambda *_: True)  # Don't actually call the collapse method
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    assert [taxa for taxa in cluster.taxa] == ['BOL', 'Bab', 'Bch', 'Bfo', 'Bfr', 'Cfu', 'Dgl', 'Edu', 'Hca', 'Hru',
                                               'Hvu', 'Lcr', 'Lla', 'Mle', 'Oma', 'Pba', 'Tin', 'Vpa']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "049088e80b31bac797a66534e518229f"
    assert cluster.taxa_sep == "-"
    assert cluster.parent is None
    assert cluster.subgroup_counter == 0
    assert cluster.cluster_score is None
    assert cluster.collapsed_genes == OrderedDict()
    assert cluster._name == "group_0"
    ids = sorted(cluster.seq_ids)
    assert ids == ['BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαE', 'BOL-PanxαF',
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
                   'Hru-PanxαC', 'Hru-PanxαD', 'Hru-PanxαE', 'Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC',
                   'Hvu-PanxβD', 'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβG', 'Hvu-PanxβH', 'Hvu-PanxβI',
                   'Hvu-PanxβJ', 'Hvu-PanxβK', 'Hvu-PanxβL', 'Hvu-PanxβM', 'Hvu-PanxβO', 'Lcr-PanxαA',
                   'Lcr-PanxαB', 'Lcr-PanxαC', 'Lcr-PanxαD', 'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαG',
                   'Lcr-PanxαH', 'Lcr-PanxαI', 'Lcr-PanxαJ', 'Lcr-PanxαK', 'Lcr-PanxαL', 'Lla-PanxαA',
                   'Lla-PanxαB', 'Lla-PanxαC', 'Mle-Panxα1', 'Mle-Panxα10A', 'Mle-Panxα11', 'Mle-Panxα12',
                   'Mle-Panxα2', 'Mle-Panxα3', 'Mle-Panxα4', 'Mle-Panxα5', 'Mle-Panxα6', 'Mle-Panxα7A',
                   'Mle-Panxα8', 'Mle-Panxα9', 'Oma-PanxαA', 'Oma-PanxαB', 'Oma-PanxαC', 'Oma-PanxαD',
                   'Pba-PanxαA', 'Pba-PanxαB', 'Pba-PanxαC', 'Pba-PanxαD', 'Pba-PanxαE', 'Pba-PanxαF',
                   'Pba-PanxαG', 'Tin-PanxαA', 'Tin-PanxαB', 'Tin-PanxαC', 'Tin-PanxαD', 'Tin-PanxαE',
                   'Tin-PanxαF', 'Vpa-PanxαA', 'Vpa-PanxαB', 'Vpa-PanxαC', 'Vpa-PanxαD', 'Vpa-PanxαE',
                   'Vpa-PanxαF', 'Vpa-PanxαG'], print(ids)
    assert cluster.seq_id_hash == "cad7ae67468eea9293c3ae2689e116ed"

    with pytest.raises(ValueError) as err:
        sim_scores = hf.get_data("cteno_sim_scores")
        sim_scores = sim_scores.ix[1:, :]
        rdmcl.Cluster(cluster.seq_ids, sim_scores)
    assert "The number of incoming sequence ids (134) does not match the expected graph size of 8911" in str(err)
    assert os.path.isfile("cad7ae67468eea9293c3ae2689e116ed.error")
    os.remove("cad7ae67468eea9293c3ae2689e116ed.error")


def test_cluster_instantiate_child(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    parent.collapsed_genes = OrderedDict([('Mle-Panxα10A', ['Mle-Panxα9'])])
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
    assert sorted(cluster.seq_ids) == child_ids
    assert cluster.seq_id_hash == hf.string2hash(", ".join(child_ids))


def test_cluster_reset_seq_ids(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    cluster.reset_seq_ids(['Dgl-PanxαE', 'Bfo-PanxαB', 'Tin-PanxαC', 'BOL-PanxαA', 'Edu-PanxαA',
                           'Bch-PanxαC', 'Vpa-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH', 'Hca-PanxαB',
                           'Oma-PanxαC', 'Bab-PanxαB', 'Mle-Panxα10A'])

    assert sorted(cluster.seq_ids) == ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                       'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A',
                                       'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    assert cluster.seq_ids_str == str(", ".join(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                                 'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A',
                                                 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB'])), print(cluster.seq_ids_str)
    assert cluster.seq_id_hash == hf.string2hash(cluster.seq_ids_str)


def test_cluster_collapse(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    assert not cluster.collapsed_genes
    cluster.collapse()
    assert cluster.collapsed_genes == OrderedDict([('Hvu-PanxβI', ['Hvu-PanxβM', 'Hvu-PanxβD', 'Hvu-PanxβC',
                                                                   'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβH',
                                                                   'Hvu-PanxβG', 'Hvu-PanxβK', 'Hvu-PanxβO',
                                                                   'Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβJ',
                                                                   'Hvu-PanxβL'])])


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
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    global_best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
    best_hits = cluster.recursive_best_hits('Bab-PanxαB', global_best_hits, ['Bab-PanxαB'])
    assert best_hits.to_csv() == """\
,psi,raw_score,score,seq1,seq2,subsmat
0,0.9841450231425388,0.9710485816574068,0.9776353490763704,Bab-PanxαB,Vpa-PanxαB,0.9748454887622982
1,1.0,0.991398404562704,1.0,Mle-Panxα9,Vpa-PanxαB,1.0
""", print(best_hits.to_csv())


def test_cluster_perturb(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    cluster.sim_scores = cluster.sim_scores[0:4]
    for i in range(0, 4):
        cluster.sim_scores.at[i, "score"] = 0.5

    assert cluster.sim_scores["score"].std() == 0

    result = cluster.perturb(cluster.sim_scores)
    assert result.iloc[0].score != 0.5
    assert round(result.iloc[0].score, 5) == 0.5


def test_cluster_get_base_cluster(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    child = rdmcl.Cluster(*hf.base_cluster_args(), parent=parent)
    grandchild = rdmcl.Cluster(*hf.base_cluster_args(), parent=child)
    assert grandchild.get_base_cluster() == parent


def test_cluster_score(hf, monkeypatch):
    monkeypatch.setattr(rdmcl.Cluster, "_score_diminishing_returns", lambda *_: 20)
    monkeypatch.setattr(rdmcl.Cluster, "_score_direct_replicate_penalty", lambda *_: 30)
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    cluster.cluster_score = 10

    assert cluster.score() == 10
    assert cluster.score(force=True) == 20
    assert cluster.score(algorithm="drp", force=True) == 30


def test_cluster_score_diminishing_returns(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())

    # No paralogs
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child.cluster_score is None

    # The second call just retrieves the attribute from the cluster saved during first call
    assert child._score_diminishing_returns() == child.cluster_score == 45.77283950617284

    # With paralogs
    child_ids = ['BOL-PanxαA', 'BOL-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child._score_diminishing_returns() == child.cluster_score == 39.713011188271594

    # Single sequence
    child_ids = ['BOL-PanxαA']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert round(child._score_diminishing_returns(), 12) == 1.847222222222

    # Include an orphan sequence
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC']
    sim_scores = hf.get_sim_scores(child_ids)
    child = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert round(child._score_diminishing_returns(), 3) == 8.575
    child.seq_ids.add("Foo-Bar3")
    assert round(child._score_diminishing_returns(), 12) == 8.510526315789

    # Edge case where child is full size of parent
    child = rdmcl.Cluster(parent.seq_ids, parent.sim_scores, parent=parent)
    assert round(child._score_diminishing_returns(), 12) == 100.877456025217


def test_get_dim_ret_base_score(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())

    # More taxa than paralogs => DRB < 0.5
    assert parent.get_dim_ret_base_score() == 0.20679012345679013

    # Num taxa == num paralogs => DRB = 0.5
    removed_taxa = ['Lla', 'Bfr', 'Oma', 'Bab', 'Bch', 'Hru', 'Cfu', 'Mle', 'Hvu', 'Lcr']
    del_list = []
    seq_ids = sorted(parent.seq_ids)
    for indx, seq_id in enumerate(seq_ids):
        if seq_id[:3] in removed_taxa:
            del_list.append(indx)
    for indx in sorted(del_list, reverse=True):
        del seq_ids[indx]
    for taxon in removed_taxa:
        del parent.taxa[taxon]

    parent.seq_ids = set(seq_ids)
    assert parent.get_dim_ret_base_score() == 0.5

    # More paralogs than taxa => DRB > 0.5
    removed_taxa = ['Tin', 'Pba', 'Vpa', 'BOL']
    del_list = []
    seq_ids = sorted(parent.seq_ids)
    for indx, seq_id in enumerate(seq_ids):
        if seq_id[:3] in removed_taxa:
            del_list.append(indx)
    for indx in sorted(del_list, reverse=True):
        del seq_ids[indx]
    for taxon in removed_taxa:
        del parent.taxa[taxon]
    parent.seq_ids = set(seq_ids)
    assert parent.get_dim_ret_base_score() == 0.7777777777777778


def test_cluster_score_direct_replicate_penalty(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())

    # No paralogs
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child.cluster_score is None

    # The second call just retrieves the attribute from the cluster saved during first call
    assert child._score_direct_replicate_penalty() == child.cluster_score == 59.475087231583096

    # With paralogs
    child_ids = ['BOL-PanxαA', 'BOL-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child._score_direct_replicate_penalty() == child.cluster_score == 42.083252786323456

    # Single sequence
    child_ids = ['BOL-PanxαA']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert round(child._score_direct_replicate_penalty(), 12) == -6.933372455078

    # Include an orphan sequence
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC']
    sim_scores = hf.get_sim_scores(child_ids)
    child = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert round(child._score_direct_replicate_penalty(), 3) == 0.309
    child.seq_ids.add("Foo-Bar3")
    assert round(child._score_direct_replicate_penalty(), 12) == 0.942727475044

    # Edge case where child is full size of parent
    child = rdmcl.Cluster(parent.seq_ids, parent.sim_scores, parent=parent)
    assert round(child._score_direct_replicate_penalty(), 12) == -263.932460129122


def test_cluster_pull_scores_subgraph(hf):
    clust = rdmcl.Cluster(*hf.base_cluster_args())
    output = str(clust.pull_scores_subgraph(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC']))
    assert output == """\
            seq1        seq2         subsmat             psi       raw_score           score
132   Hru-PanxαA  Lcr-PanxαH  0.974895831252  0.952850560876  0.965773924299  0.968282250139
165   Hru-PanxαA  Tin-PanxαC  0.968036771413  0.937704276879  0.958388562706  0.958937023053
6851  Lcr-PanxαH  Tin-PanxαC  0.984286835987  0.968945697684  0.975097822535  0.979684494496""", print(output)


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
    sim_scores.at[indx.index[0], "raw_score"] = 0.7

    indx = sim_scores[(sim_scores.seq1 == 'Mle-Panxα10A') & (sim_scores.seq2 == 'Vpa-PanxαB')]
    sim_scores.at[indx.index[0], "raw_score"] = 0.8

    indx = sim_scores[(sim_scores.seq2 == 'Vpa-PanxαB') & (sim_scores.seq1 == 'BOL-PanxαA')]
    sim_scores.at[indx.index[0], "raw_score"] = 0.99

    indx = sim_scores[(sim_scores.seq2 == 'Mle-Panxα9') & (sim_scores.seq1 == 'BOL-PanxαA')]
    sim_scores.at[indx.index[0], "raw_score"] = 0.99

    indx = sim_scores[(sim_scores.seq1 == 'Mle-Panxα9') & (sim_scores.seq2 == 'Vpa-PanxαB')]
    sim_scores.at[indx.index[0], "raw_score"] = 0.97

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
    assert "TERMINATED: Spinning off cliques would orphan Bch-PanxαB and Lcr-PanxαI" in log_file.read(), print(log_file.read())


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
    assert sorted(cliques[0].seq_ids) == ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4']
    assert sorted(cliques[1].seq_ids) == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαA']

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
    assert sorted(cliques[0].seq_ids) == ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4', 'Vpa-PanxαA']
    assert sorted(cliques[1].seq_ids) == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
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
    assert sorted(new_clusters[0].seq_ids) == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    assert sorted(new_clusters[1].seq_ids) == ['Edu-PanxαG', 'Lcr-PanxαA', 'Pba-PanxαE']
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


def test_cluster2database(hf, monkeypatch):
    monkeypatch.setattr(rdmcl.Cluster, "score", lambda *_: 20)
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
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
    assert response[0][0] == "cad7ae67468eea9293c3ae2689e116ed"                 # hash
    assert 'BOL-PanxαA, BOL-PanxαB, BOL-PanxαC, BOL-PanxαD' in response[0][1]   # seq_ids
    assert response[0][2] == '>Seq1\nMPQQCS-SS\n>Seq2\nMPQICMAAS'               # alignment
    assert 'Hca-PanxαG,Lla-PanxαC,0.291912106162889' in response[0][3]          # graph
    assert response[0][4] == '20'                                               # score
    connect.close()


# #########  PSI-PRED  ########## #
def test_mc_psi_pred(hf, monkeypatch):
    outdir = br.TempDir()
    with open(os.path.join(hf.resource_path, "psi_pred", "BOL-PanxαB.ss2"), "r") as ofile:
        ss2_file = ofile.read()
    monkeypatch.setattr(rdmcl, "run_psi_pred", lambda *_: ss2_file)
    seq_rec = hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"]
    rdmcl.mc_psi_pred(seq_rec, [outdir.path])
    with open(os.path.join(outdir.path, "BOL-PanxαB.ss2"), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open(os.path.join(outdir.path, "BOL-PanxαB.ss2"), "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.mc_psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [outdir.path])
    with open(os.path.join(outdir.path, "BOL-PanxαB.ss2"), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "af9666d37426caa2bbf6b9075ce8df96", print(output)


def test_run_psi_pred(hf, monkeypatch, capsys):
    tempdir = br.TempDir()
    monkeypatch.setattr(br, "TempDir", lambda *_: tempdir)

    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Bab-PanxαA")
    ss2_file = rdmcl.run_psi_pred(seqbuddy.records[0])
    assert hf.string2hash(ss2_file) == "b50f39dc22e4d16be325efdd14f7900d"

    monkeypatch.setattr(shutil, "which", lambda *_: False)
    monkeypatch.setattr(rdmcl, "Popen", MockPopen)
    rdmcl.run_psi_pred(seqbuddy.records[0])
    out, err = capsys.readouterr()
    assert "rdmcl/psipred/bin/seq2mtx sequence.fa" in out
    assert "rdmcl/psipred/bin/psipred" in out
    assert "rdmcl/psipred/bin/psipass2" in out


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
    ss2_1 = os.path.join(hf.resource_path, "psi_pred", "Mle-Panxα10A.ss2")
    ss2_1 = pd.read_csv(ss2_1, comment="#", header=None, delim_whitespace=True)
    ss2_1.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]

    ss2_2 = os.path.join(hf.resource_path, "psi_pred", "Mle-Panxα8.ss2")
    ss2_2 = pd.read_csv(ss2_2, comment="#", header=None, delim_whitespace=True)
    ss2_2.columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]

    assert rdmcl.compare_psi_pred(ss2_1, ss2_2) == 0.691672882672883


# #########  Orthogroup caller  ########## #
def test_orthogroup_caller(hf):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    seq_ids = ['BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαF', 'BOL-PanxαG', 'BOL-PanxαH', 'Lcr-PanxαA',
               'Lcr-PanxαB', 'Lcr-PanxαD', 'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαH', 'Lcr-PanxαI',
               'Lcr-PanxαK', 'Mle-Panxα1', 'Mle-Panxα2', 'Mle-Panxα4', 'Mle-Panxα5', 'Mle-Panxα6',
               'Mle-Panxα7A', 'Mle-Panxα8', 'Mle-Panxα10A', 'Vpa-PanxαA',
               'Vpa-PanxαC', 'Vpa-PanxαF']

    seqbuddy = rdmcl.Sb.SeqBuddy(hf.resource_path + "Cteno_pannexins.fa")
    rdmcl.Sb.pull_recs(seqbuddy, "^%s$" % "$|^".join(seq_ids))

    # psi_pred_ss2_paths = hf.get_data("ss2_paths")
    # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, psi_pred_ss2_paths, broker)

    cluster = rdmcl.Cluster(seq_ids, hf.get_db_graph("94650a087fe83a6ad6c6584e76025b68", broker))
    cluster.collapsed_genes = OrderedDict([('Lcr-PanxαA', ['Lcr-PanxαL']), ('Lcr-PanxαE', ['Lcr-PanxαJ']),
                                           ('Mle-Panxα10A', ['Mle-Panxα9']),
                                           ('Vpa-PanxαC', ['Vpa-PanxαG', 'Vpa-PanxαB']),
                                           ('Mle-Panxα6', ['Mle-Panxα12']), ('BOL-PanxαC', ['BOL-PanxαA'])])
    cluster_list = []
    outdir = br.TempDir()
    outdir.subdir("progress")
    outdir.subdir("alignments")
    outdir.subdir("mcmcmc")
    outdir.subdir("psi_pred")
    outdir.subdir("sim_scores")
    outdir.subdir("hmm")

    progress = rdmcl.Progress(os.path.join(outdir.path, "progress"), cluster)
    psi_pred_ss2_dfs = OrderedDict()
    for rec in seq_ids:
        psi_pred_ss2_dfs[rec] = os.path.join(hf.resource_path, "psi_pred", "%s.ss2" % rec)

    steps = 10
    r_seed = 2
    orthogroups = rdmcl.orthogroup_caller(cluster, cluster_list, seqbuddy, broker, progress, outdir.path,
                                          psi_pred_ss2_dfs, steps=steps, r_seed=r_seed)

    expected = [['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4', 'Vpa-PanxαA'],
                ['Lcr-PanxαA', 'Mle-Panxα5', 'Vpa-PanxαF'],
                ['BOL-PanxαD', 'Lcr-PanxαD', 'Mle-Panxα2'],
                ['BOL-PanxαC', 'Lcr-PanxαE', 'Mle-Panxα6'],
                ['Lcr-PanxαK', 'Mle-Panxα7A'],
                ['Lcr-PanxαH', 'Mle-Panxα10A'],
                ['Lcr-PanxαB', 'Mle-Panxα1'],
                ['BOL-PanxαH', 'Mle-Panxα8'],
                ['BOL-PanxαG', 'Lcr-PanxαF'],
                ['Vpa-PanxαC']]

    orthogroup_seqs = [sorted(clust.seq_ids) for clust in orthogroups]
    assert sorted(cluster.seq_ids) in orthogroup_seqs
    for clust in expected:
        assert clust in orthogroup_seqs, print(orthogroup_seqs)


# #########  Miscellaneous  ########## #
def test_progress(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    tmpdir = br.TempDir()
    progress = rdmcl.Progress(tmpdir.path, cluster)
    assert progress.outdir == tmpdir.path
    assert os.path.isfile("{0}{1}.progress".format(tmpdir.path, hf.sep))
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # The dictionary is not static, so just sort the string:
        # {"placed": 0, "mcl_runs": 0, "total": 119}
        assert "".join(sorted(ifile.read())) == '     """""",,00134:::_aaccdelllmnoprsttu{}'

    progress.update("mcl_runs", 2)
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # {"placed": 0, "mcl_runs": 2, "total": 119}
        assert "".join(sorted(ifile.read())) == '     """""",,01234:::_aaccdelllmnoprsttu{}'

    json = progress.read()
    assert json["mcl_runs"] == 2
    assert json["placed"] == 0
    assert json["total"] == 134

    assert str(progress) == "MCL runs processed: 2. Sequences placed: 0/134. Run time: "


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


def test_heartbeat_init():
    tmpdir = br.TempDir()
    hbdb_path = tmpdir.subfile("hbdb.sqlite")
    heartbeat = rdmcl.HeartBeat(hbdb_path, 60)
    assert heartbeat.hbdb_path == hbdb_path
    assert heartbeat.pulse_rate == 60
    assert heartbeat.id is None
    assert heartbeat.running_process is None
    assert heartbeat.thread_type == "master"
    assert not heartbeat.dummy

    conn = sqlite3.connect(hbdb_path)
    cursor = conn.cursor()
    tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    assert tables == [("heartbeat",), ('sqlite_sequence',)]

    # Instantiate a second object to ensure 'table exists' error is handled
    rdmcl.HeartBeat(hbdb_path, 60)
    tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
    assert tables == [("heartbeat",), ('sqlite_sequence',)]


def test_heartbeat_dummy():
    tmpdir = br.TempDir()
    hbdb_path = tmpdir.subfile("hbdb.sqlite")
    heartbeat = rdmcl.HeartBeat(hbdb_path, 60, dummy=True)

    conn = sqlite3.connect(hbdb_path)
    cursor = conn.cursor()
    assert not cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()

    heartbeat._run("foo")
    assert heartbeat.id is None
    assert heartbeat.running_process is None
    assert not cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()

    heartbeat.start()
    assert heartbeat.id is None
    assert heartbeat.running_process is None
    assert not cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()

    heartbeat.end()
    assert heartbeat.id is None
    assert heartbeat.running_process is None
    assert not cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()


def test_heartbeat_run(monkeypatch):
    tmpdir = br.TempDir()
    checkfile = tmpdir.subfile("checkfile")
    with open(checkfile, "w") as ofile:
        ofile.write("Running")

    hbdb_path = tmpdir.subfile("hbdb.sqlite")
    heartbeat = rdmcl.HeartBeat(hbdb_path, 0)
    heartbeat.id = 1
    conn = sqlite3.connect(hbdb_path)
    cursor = conn.cursor()
    start_time = round(rdmcl.time.time()) - 1
    cursor.execute("INSERT INTO heartbeat (thread_type, pulse) " 
                   "VALUES (?, ?)", ("master", start_time,))
    conn.commit()

    monkeypatch.setattr(rdmcl.time, "sleep", mock_keyboardinterupt)
    heartbeat._run(checkfile)

    query = cursor.execute("SELECT * FROM heartbeat").fetchall()
    assert len(query) == 1
    assert query[0][2] > start_time

    monkeypatch.setattr(rdmcl.time, "sleep", lambda *_: open(checkfile, "w").close())
    heartbeat._run(checkfile)
    query2 = cursor.execute("SELECT * FROM heartbeat").fetchall()
    assert query[0][2] <= query2[0][2]


def test_heartbeat_start(monkeypatch):
    tmpdir = br.TempDir()
    hbdb_path = tmpdir.subfile("hbdb.sqlite")
    heartbeat = rdmcl.HeartBeat(hbdb_path, 0, thread_type="worker")

    conn = sqlite3.connect(hbdb_path)
    cursor = conn.cursor()
    monkeypatch.setattr(rdmcl.time, "sleep", mock_keyboardinterupt)
    monkeypatch.setattr(rdmcl.time, "time", lambda *_: 123456)

    heartbeat.start()
    assert heartbeat.check_file.read() == "Running"  # This is a race, but should be fine.
    while heartbeat.running_process.is_alive():
        pass
    assert heartbeat.check_file.read() == ""
    assert heartbeat.id == 1
    query = cursor.execute("SELECT * FROM heartbeat").fetchall()
    assert len(query) == 1
    assert query[0][0] == 1
    assert query[0][1] == "worker"
    assert query[0][2] == 123456

    monkeypatch.setattr(rdmcl.HeartBeat, "end", lambda *_: True)
    heartbeat.start()
    assert heartbeat.id == 2
    query = cursor.execute("SELECT * FROM heartbeat").fetchall()
    assert len(query) == 2
    assert query[1][0] == 2
    assert query[1][1] == "worker"
    assert query[1][2] == 123456


def test_heartbeat_end():
    tmpdir = br.TempDir()
    hbdb_path = tmpdir.subfile("hbdb.sqlite")
    heartbeat = rdmcl.HeartBeat(hbdb_path, 0, thread_type="master")

    heartbeat.id = 1
    conn = sqlite3.connect(hbdb_path)
    cursor = conn.cursor()
    cursor.execute("INSERT INTO heartbeat (thread_type, pulse) " 
                   "VALUES (?, ?)", ("master", 123456,))
    conn.commit()

    # Not running yet
    heartbeat.end()
    assert cursor.execute("SELECT thread_type FROM heartbeat").fetchone() == ("master",)

    class MockProcess(object):
        def __init__(self):
            self.alive = True

        def is_alive(self):
            if self.alive:
                self.alive = False
                return True
            else:
                self.alive = True
                return False

    heartbeat.running_process = MockProcess()
    heartbeat.end()
    assert heartbeat.running_process is None
    assert not cursor.execute("SELECT thread_type FROM heartbeat").fetchone()
    assert heartbeat.id is None


# ################ SCORING FUNCTIONS ################ #
def test_mc_score_sequences1(hf):
    outfile = br.TempFile()
    alb_obj = hf.get_data("cteno_panxs_aln")
    # Grab subset of alignment for manual calculation
    # Bfo-PanxαA   SQMWSQ--DDA
    # Bfr-PanxαD   V--RQIVVGGP
    alb_obj = rdmcl.Alb.extract_regions(alb_obj, "105:115")

    ss2_dfs = hf.get_data("ss2_dfs")
    ss2_dfs = {"Bfo-PanxαA": ss2_dfs["Bfo-PanxαA"], "Bfr-PanxαD": ss2_dfs["Bfr-PanxαD"]}
    # Update the ss2 dfs according to the alignment subsequences extracted
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
47     1  S  H      0.034       0.966       0.003
48     2  Q  H      0.071       0.926       0.004
49     3  M  H      0.371       0.649       0.003
50     4  W  C      0.802       0.211       0.004
51     5  S  C      0.852       0.151       0.010
52     6  Q  C      0.765       0.253       0.009
53     9  D  C      0.733       0.283       0.011
54    10  D  C      0.890       0.126       0.014
55    11  A  C      0.914       0.085       0.028
    '''
    '''
    indx aa ss  coil_prob  helix_prob  sheet_prob
44     1  V  E      0.136       0.196       0.556
45     4  R  E      0.178       0.357       0.553
46     5  Q  H      0.157       0.530       0.525
47     6  I  H      0.114       0.771       0.319
48     7  V  H      0.217       0.675       0.212
49     8  V  C      0.461       0.443       0.166
50     9  G  C      0.837       0.040       0.077
51    10  G  C      0.940       0.008       0.063
52    11  P  C      0.606       0.015       0.402
    '''

    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].iloc[47:56]
    for indx, new in [(47, 1), (48, 2), (49, 3), (50, 4), (51, 5), (52, 6), (53, 9), (54, 10), (55, 11)]:
        ss2_dfs["Bfo-PanxαA"].at[indx, "indx"] = new

    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].iloc[44:53]
    for indx, new in [(44, 1), (45, 4), (46, 5), (47, 6), (48, 7), (49, 8), (50, 9), (51, 10), (52, 11)]:
        ss2_dfs["Bfr-PanxαD"].at[indx, "indx"] = new

    gap_open = -5
    gap_extend = 0

    # For score, subsmat = -0.363
    rdmcl.mc_score_sequences([("Bfo-PanxαA", "Bfr-PanxαD", ss2_dfs["Bfo-PanxαA"], ss2_dfs["Bfr-PanxαD"])],
                             [alb_obj, gap_open, gap_extend, outfile.path])

    assert outfile.read() == "\nBfo-PanxαA,Bfr-PanxαD,-0.3627272727272728,0.4183636363636363"


def test_mc_score_sequences2(monkeypatch):
    results_file = br.TempFile()
    seq_pairs = [("seq1", "seq2", "fakedf1", "fakedf2"), ("seq1", "seq3",  "fakedf1", "fakedf3"),
                 ("seq2", "seq3",  "fakedf2", "fakedf3")]
    alb_obj = rdmcl.Alb.AlignBuddy("""\
>seq1
MP-QMSASWI
>Seq2
MPPQISAS-I
>Seq3
MP-QISGAWI
""")

    args = [alb_obj, -5, 0, results_file.path]

    monkeypatch.setattr(rdmcl, "compare_pairwise_alignment", lambda *_: "subs_mat_score")
    monkeypatch.setattr(rdmcl, "compare_psi_pred", lambda *_: "ss_score")

    rdmcl.mc_score_sequences(seq_pairs, args)

    results = results_file.read()
    assert results == """
seq1,seq2,subs_mat_score,ss_score
seq1,seq3,subs_mat_score,ss_score
seq2,seq3,subs_mat_score,ss_score""", print(results)


def test_compare_pairwise_alignment():
    alignbuddy = rdmcl.Alb.AlignBuddy("""\
>seq1
MP--QMSASWI
>Seq2
MPPIQISAS-I
""")
    subs_mat_score = rdmcl.compare_pairwise_alignment(alignbuddy, -5, -1)
    assert subs_mat_score == 0.40982529375386517


def test_mc_create_all_by_all_scores(capsys, monkeypatch):
    monkeypatch.setattr(rdmcl, "retrieve_all_by_all_scores", lambda *args, **kwargs: print(args, kwargs))
    rdmcl.mc_create_all_by_all_scores("seqbuddy", ["arg1", "arg2"])
    out, err = capsys.readouterr()
    assert out == "('seqbuddy', 'arg1', 'arg2') {'quiet': True}\n"


def test_retrieve_all_by_all_scores_single(hf):
    sql_broker = helpers.SQLiteBroker(os.path.join(hf.resource_path, "db.sqlite"))
    sql_broker.start_broker()

    seqbuddy = rdmcl.Sb.pull_recs(hf.get_data("cteno_panxs"), "Mle-Panxα5")
    sim_scores, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, "psi_pred_files", sql_broker)
    assert len(alignbuddy.records()) == 1
    assert sim_scores.empty
    assert list(sim_scores.columns) == ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
    sql_broker.close()


def test_retrieve_all_by_all_scores_from_db(hf, monkeypatch):
    sql_broker = helpers.SQLiteBroker(os.path.join(hf.resource_path, "db.sqlite"))
    sql_broker.start_broker()

    seqbuddy = rdmcl.Sb.pull_recs(hf.get_data("cteno_panxs"), "Bfo-PanxαF|Hca-PanxαD|Mle-Panxα6")

    # These patches shouldn't execute, they're here to make sure the control flow doesn't go passed the query statement
    monkeypatch.setattr(rdmcl, "WorkerJob", mock_keyboardinterupt)
    monkeypatch.setattr(rdmcl.AllByAllScores, "create", mock_keyboardinterupt)
    # sim_scores, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, hf.get_data("ss2_paths"), sql_broker)
    sim_scores, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, "psi_pred_files", sql_broker)
    assert len(alignbuddy.records()) == 3
    assert str(sim_scores) == """\
         seq1        seq2         subsmat             psi       raw_score           score
0  Bfo-PanxαF  Hca-PanxαD  0.117138287495  0.000000000000  0.820955423002  0.081996801246
1  Bfo-PanxαF  Mle-Panxα6  1.000000000000  1.000000000000  0.913305173757  1.000000000000
2  Hca-PanxαD  Mle-Panxα6  0.000000000000  0.192926181927  0.814352991900  0.057877854578""", print(sim_scores)
    sql_broker.close()


def test_retrieve_all_by_all_scores_feed_worker(hf, monkeypatch):
    sql_broker = helpers.SQLiteBroker(os.path.join(hf.resource_path, "db.sqlite"))
    sql_broker.start_broker()

    seqbuddy = rdmcl.Sb.pull_recs(hf.get_data("cteno_panxs"), "Mle")
    tmpfile = br.TempFile()

    class MockWorkerJob(object):
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

        @staticmethod
        def run():
            return "worker_sim_scores", "worker_alignment"

    monkeypatch.setattr(rdmcl, "WorkerJob", MockWorkerJob)
    monkeypatch.setattr(rdmcl, "WORKER_DB", tmpfile.path)
    monkeypatch.setattr(rdmcl, "MIN_SIZE_TO_WORKER", 1)

    sim_scores, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, "psi_pred_files", sql_broker)
    assert sim_scores == "worker_sim_scores"
    assert alignbuddy == "worker_alignment"


def test_retrieve_all_by_all_scores_new_run(hf, monkeypatch):
    sql_broker = helpers.SQLiteBroker(os.path.join(hf.resource_path, "db.sqlite"))
    sql_broker.start_broker()

    monkeypatch.setattr(rdmcl, "MIN_SIZE_TO_WORKER", 1000)
    monkeypatch.setattr(rdmcl.AllByAllScores, "create", lambda *_, **__: ["create_sim_scores", "create_alignment"])
    seqbuddy = rdmcl.Sb.pull_recs(hf.get_data("cteno_panxs"), "Mle")

    sim_scores, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, "psi_pred_files", sql_broker)
    assert sim_scores == "create_sim_scores"
    assert alignbuddy == "create_alignment"
    sql_broker.close()


def test_allbyallscores_init(hf):
    seqbuddy = rdmcl.Sb.pull_recs(hf.get_data("cteno_panxs"), "Mle")
    all_by_all = rdmcl.AllByAllScores(seqbuddy, "psipred", "sqlbroker")
    assert type(all_by_all.seqbuddy) == rdmcl.Sb.SeqBuddy
    assert all_by_all.seqbuddy != seqbuddy
    assert str(all_by_all.seqbuddy.records) == str(seqbuddy.records)
    for rec in all_by_all.seqbuddy.records:
        assert rec.id.startswith("Mle")
    assert all_by_all.psi_pred_ss2 == "psipred"
    assert all_by_all.sql_broker == "sqlbroker"
    assert not all_by_all.quiet


def test_allbyallscores_create(hf):
    tmpdir = br.TempDir()
    sql_broker = helpers.SQLiteBroker(os.path.join(tmpdir.path, "db.sqlite"))
    sql_broker.create_table("data_table", ["hash TEXT PRIMARY KEY", "seq_ids TEXT", "alignment TEXT",
                                           "graph TEXT", "cluster_score TEXT"])
    sql_broker.start_broker()

    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    psi_pred_files = [(rec.id, os.path.join(hf.resource_path, "psi_pred", "%s.ss2" % rec.id))
                      for rec in seqbuddy.records]
    psi_pred_files = OrderedDict(psi_pred_files)

    all_by_all_obj = rdmcl.AllByAllScores(seqbuddy, psi_pred_files, sql_broker)
    sim_scores, alignbuddy = all_by_all_obj.create()
    assert len(sim_scores.index) == 66  # This is for 12 starting sequences --> (a * (a - 1)) / 2
    compare = sim_scores.loc[:][(sim_scores['seq1'] == "Mle-Panxα2") & (sim_scores['seq2'] == "Mle-Panxα12")]
    assert "Mle-Panxα2  Mle-Panxα12  0.332679039959  0.37027112509  0.53661224001" in str(compare), print(compare)
    assert len(alignbuddy.records()) == 12
    sql_broker.close()


def test_worker_update_psipred(hf):
    align = hf.get_data("cteno_panxs_aln")
    align = rdmcl.Alb.extract_regions(align, "105:115")
    align = rdmcl.Alb.pull_records(align, "Bfo-PanxαA|Bfr-PanxαD")
    ss2_dfs = hf.get_data("ss2_dfs")
    ss2_dfs = {"Bfo-PanxαA": ss2_dfs["Bfo-PanxαA"].iloc[47:56], "Bfr-PanxαD": ss2_dfs["Bfr-PanxαD"].iloc[44:53]}

    for indx, new in [(47, 1), (48, 2), (49, 3), (50, 4), (51, 5), (52, 6), (53, 7), (54, 8), (55, 9)]:
        ss2_dfs["Bfo-PanxαA"].at[indx, "indx"] = new
    ss2_dfs["Bfo-PanxαA"] = ss2_dfs["Bfo-PanxαA"].reset_index(drop=True)

    for indx, new in [(44, 1), (45, 2), (46, 3), (47, 4), (48, 5), (49, 6), (50, 7), (51, 8), (52, 9)]:
        ss2_dfs["Bfr-PanxαD"].at[indx, "indx"] = new
    ss2_dfs["Bfr-PanxαD"] = ss2_dfs["Bfr-PanxαD"].reset_index(drop=True)

    ss2_dfs = rdmcl.update_psipred(align, ss2_dfs, "msa")
    assert str(ss2_dfs["Bfo-PanxαA"]) == """\
   indx aa ss  coil_prob  helix_prob  sheet_prob
0     0  S  H      0.034       0.966       0.003
1     1  Q  H      0.071       0.926       0.004
2     2  M  H      0.371       0.649       0.003
3     3  W  C      0.802       0.211       0.004
4     4  S  C      0.852       0.151       0.010
5     5  Q  C      0.765       0.253       0.009
6     8  D  C      0.733       0.283       0.011
7     9  D  C      0.890       0.126       0.014
8    10  A  C      0.914       0.085       0.028"""

    align = rdmcl.Alb.trimal(align, "all")
    ss2_dfs = rdmcl.update_psipred(align, ss2_dfs, "trimal")
    assert str(ss2_dfs["Bfo-PanxαA"]) == """\
   indx aa ss  coil_prob  helix_prob  sheet_prob
0     0  S  H      0.034       0.966       0.003
1     3  W  C      0.802       0.211       0.004
2     4  S  C      0.852       0.151       0.010
3     5  Q  C      0.765       0.253       0.009
4     8  D  C      0.733       0.283       0.011
5     9  D  C      0.890       0.126       0.014
6    10  A  C      0.914       0.085       0.028""", print(str(ss2_dfs["Bfo-PanxαA"]))

    with pytest.raises(ValueError) as err:
        rdmcl.update_psipred(align, ss2_dfs, "foo")

    assert "Unrecognized mode 'foo': select from ['msa', 'trimal']" in str(err)


def test_trimal():
    align = rdmcl.Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STP---YWAILP
""", in_format="fasta")

    seqbuddy = rdmcl.Sb.SeqBuddy(rdmcl.Alb.make_copy(align).records(), in_format="fasta")
    seqbuddy = rdmcl.Sb.clean_seq(seqbuddy)

    # Don't modify if any sequence is reduced to nothing
    trimal = rdmcl.trimal(seqbuddy, [0.3], rdmcl.Alb.make_copy(align))
    assert str(trimal) == str(align)

    align = rdmcl.Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STPTC-YWAILP
""", in_format="fasta")

    # Don't modify if average sequence length is reduced by more than half
    trimal = rdmcl.trimal(seqbuddy, [0.3], rdmcl.Alb.make_copy(align))
    assert str(trimal) == str(align)

    # Remove some gaps
    trimal = rdmcl.trimal(seqbuddy, ["all", 0.3, 0.55, "clean"], rdmcl.Alb.make_copy(align))
    assert str(trimal) == """\
>A
MSTGTC----
>B
M---TC----
>C
M---TCAILP
>D
-STPTCAILP
"""


def test_prepare_all_by_all(hf):
    cpus = 24
    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy = rdmcl.Sb.pull_recs(seqbuddy, "Oma")  # Only 4 records, which means 6 comparisons
    ss2_dfs = hf.get_data("ss2_dfs")

    data_len, data = rdmcl.prepare_all_by_all(seqbuddy, ss2_dfs, cpus)
    assert data_len == 6
    assert len(data) == 6
    assert data[0][0][0:2] == ('Oma-PanxαA', 'Oma-PanxαB')

    seqbuddy = hf.get_data("cteno_panxs")  # 134 records = 8911 comparisons
    data_len, data = rdmcl.prepare_all_by_all(seqbuddy, ss2_dfs, cpus)
    assert data_len == 8911
    assert len(data[0]) == int(ceil(data_len / cpus))


def test_set_final_sim_scores(hf):
    sim_scores = hf.get_data("cteno_sim_scores")
    sim_scores = sim_scores.iloc[:10]
    sim_scores = rdmcl.set_final_sim_scores(sim_scores)
    assert str(sim_scores) == """\
         seq1         seq2         subsmat             psi       raw_score           score
0  Hca-PanxαG   Lla-PanxαC  0.239664985272  0.000000000000  0.320014096667  0.167765489691
1  Hca-PanxαG   Mle-Panxα1  0.000000000000  0.453888468075  0.301324591814  0.136166540423
2  Hca-PanxαG   Mle-Panxα2  0.129888458521  0.573451478794  0.350239526820  0.262957364603
3  Hca-PanxαG   Mle-Panxα3  0.374897044166  0.635128529572  0.425122352909  0.452966489787
4  Hca-PanxαG   Mle-Panxα4  1.000000000000  1.000000000000  0.638192277715  1.000000000000
5  Hca-PanxαG   Mle-Panxα5  0.189100677982  0.556072403045  0.364911258772  0.299192195501
6  Hca-PanxαG   Mle-Panxα6  0.424713073878  0.597077241840  0.434979887508  0.476422324267
7  Hca-PanxαG  Mle-Panxα7A  0.190522025243  0.681611915423  0.378627780350  0.337848992297
8  Hca-PanxαG   Mle-Panxα8  0.290356199296  0.511684844001  0.388444197191  0.356754792708
9  Hca-PanxαG   Mle-Panxα9  0.483807414359  0.580623914006  0.449716964507  0.512852364253""", print(sim_scores)


def test_workerjob_init(hf, monkeypatch):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    seqbuddy = hf.get_data("cteno_panxs")
    seq_ids = hf.get_data("cteno_ids")
    seq_ids_hash = hf.string2hash(", ".join(seq_ids))
    sql_broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    sql_broker.start_broker()
    job_id = "".join([str(x) for x in [seq_ids_hash, rdmcl.GAP_OPEN, rdmcl.GAP_EXTEND,
                                       rdmcl.ALIGNMETHOD, rdmcl.ALIGNPARAMS, rdmcl.TRIMAL]])
    job_id = hf.string2hash(job_id)

    worker = rdmcl.WorkerJob(seqbuddy, sql_broker)
    assert worker.seqbuddy == seqbuddy
    assert worker.seq_ids == seq_ids
    assert worker.seq_id_hash == seq_ids_hash
    assert worker.job_id == job_id
    assert worker.sql_broker == sql_broker
    assert worker.heartbeat.hbdb_path == rdmcl.HEARTBEAT_DB
    assert worker.heartbeat.pulse_rate == rdmcl.MASTER_PULSE


def test_workerjob_run(hf, monkeypatch, capsys):
    temp_dir = br.TempDir()
    hb_db = temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)

    hb_con = sqlite3.connect(os.path.join(temp_dir.path, "heartbeat_db.sqlite"))
    hb_cursor = hb_con.cursor()

    hb_cursor.execute("DELETE FROM heartbeat")
    hb_con.commit()

    seqbuddy = hf.get_data("cteno_panxs")

    rdmcl.HEARTBEAT_DB = hb_db
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    monkeypatch.setattr(rdmcl.WorkerJob, "queue_job", lambda *_: print("queue_job()"))

    worker = rdmcl.WorkerJob(seqbuddy, "sql_broker")
    # No workers available
    assert not worker.run()
    out, err = capsys.readouterr()
    assert "queue_job" not in out

    # Return DB result
    monkeypatch.setattr(rdmcl.WorkerJob, "pull_from_db", lambda *_: "db_result")
    hb_cursor.execute("INSERT INTO heartbeat (thread_type, pulse) VALUES ('worker', %s)" % 10e12)
    hb_con.commit()
    assert worker.run() == "db_result"
    out, err = capsys.readouterr()
    assert "queue_job" in out

    # Job finished
    def mock_pull_from_db(self):
        self.running = False
        return False

    monkeypatch.setattr(rdmcl.WorkerJob, "pull_from_db", mock_pull_from_db)
    monkeypatch.setattr(rdmcl.WorkerJob, "check_finished", lambda *_: True)
    monkeypatch.setattr(rdmcl.WorkerJob, "process_finished", lambda *_: "process_finished")

    assert worker.run() == "process_finished"

    # Job died
    monkeypatch.setattr(rdmcl.WorkerJob, "check_finished", lambda *_: False)
    monkeypatch.setattr(rdmcl, "random", lambda *_: 0.99)
    monkeypatch.setattr(rdmcl.WorkerJob, "check_if_active", lambda *_: False)
    assert not worker.run()

    # Go to the end and hit pause
    monkeypatch.setattr(rdmcl, "random", lambda *_: 0.01)
    assert not worker.run()


def test_workerjob_queue_job(hf, monkeypatch):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    temp_dir = br.TempDir()
    work_db = temp_dir.copy_to(os.path.join(hf.resource_path, "work_db.sqlite"))
    rdmcl.WORKER_DB = work_db
    rdmcl.WORKER_OUT = temp_dir.path

    work_con = sqlite3.connect(work_db)
    work_cursor = work_con.cursor()

    seqbuddy = hf.get_data("cteno_panxs")
    worker = rdmcl.WorkerJob(seqbuddy, "sql_broker")

    assert not worker.queue_job()
    assert os.path.isfile(os.path.join(temp_dir.path, "%s.seqs" % worker.job_id))
    queue = work_cursor.execute("SELECT * FROM queue").fetchall()
    assert queue == [('a2aaca4f79bd56fbf8debfdc281660fd', '', 'clustalo', '',
                      'gappyout 0.5 0.75 0.9 0.95 clean', -5.0, 0.0)], print(queue)


def test_workerjob_pull_from_db(hf, monkeypatch):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    temp_dir = br.TempDir()
    broker_db = temp_dir.copy_to(os.path.join(hf.resource_path, "db.sqlite"))
    sql_broker = helpers.SQLiteBroker(broker_db)
    sql_broker.start_broker()

    work_db = temp_dir.copy_to(os.path.join(hf.resource_path, "work_db.sqlite"))
    rdmcl.WORKER_DB = work_db
    rdmcl.WORKER_OUT = temp_dir.path

    work_con = sqlite3.connect(work_db)
    work_cursor = work_con.cursor()

    seqbuddy = hf.get_data("cteno_panxs")
    worker = rdmcl.WorkerJob(seqbuddy, sql_broker)
    worker.heartbeat.id = 1

    # No record in DB
    assert not worker.pull_from_db()

    # Get record
    graph = """\
Hca-PanxαG,Lla-PanxαC,0.239665,0.000000,0.320014,0.167765
Hca-PanxαG,Mle-Panxα1,0.000000,0.453888,0.301325,0.136167
"""
    alignment = """\
>Hca-PanxαG
MTGLILIL
>Lla-PanxαC
MTGLILIL
>Mle-Panxα1
MTGLILIL
"""
    sql_broker.query("INSERT INTO data_table (hash, graph, alignment) VALUES (?, ?, ?)", (worker.seq_id_hash,
                                                                                          graph, alignment,))
    work_cursor.execute("INSERT INTO complete (hash) VALUES (?)", (worker.job_id,))
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", (worker.job_id, 1))
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", (worker.job_id, 2))
    work_con.commit()

    pulled_data = worker.pull_from_db()
    assert str(pulled_data[0]) == """\
         seq1        seq2   subsmat       psi  raw_score     score
0  Hca-PanxαG  Lla-PanxαC  0.239665  0.000000   0.320014  0.167765
1  Hca-PanxαG  Mle-Panxα1  0.000000  0.453888   0.301325  0.136167"""
    assert str(pulled_data[1]) == alignment

    assert len(work_cursor.execute("SELECT * FROM waiting").fetchall()) == 1

    worker.heartbeat.id = 2
    worker.pull_from_db()
    assert not work_cursor.execute("SELECT * FROM waiting").fetchall()


def test_workerjob_check_finished(hf, monkeypatch):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    temp_dir = br.TempDir()
    work_db = temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    rdmcl.WORKER_DB = work_db
    rdmcl.WORKER_OUT = temp_dir.path

    work_con = sqlite3.connect(work_db)
    work_cursor = work_con.cursor()

    seqbuddy = hf.get_data("cteno_panxs")
    worker = rdmcl.WorkerJob(seqbuddy, "sql_broker")
    worker.heartbeat.id = 1

    work_cursor.execute("INSERT INTO queue (hash, psi_pred_dir, align_m, align_p, trimal, "
                        "gap_open, gap_extend) VALUES (?, ?, ?, ?, ?, ?, ?)",
                        ('a2aaca4f79bd56fbf8debfdc281660fd', '', 'clustalo', '',
                         'gappyout 0.5 0.75 0.9 0.95 clean', -5.0, 0.0,))
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", ('a2aaca4f79bd56fbf8debfdc281660fd', 1,))
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", ('a2aaca4f79bd56fbf8debfdc281660fd', 2,))
    work_con.commit()

    # Not finished
    assert not worker.check_finished()
    assert work_cursor.execute("SELECT  COUNT(*) FROM waiting").fetchone()[0] == 2

    # Other jobs waiting
    work_cursor.execute("INSERT INTO complete (hash) "
                        "VALUES (?)", ('a2aaca4f79bd56fbf8debfdc281660fd',))
    work_con.commit()
    assert worker.check_finished()
    assert work_cursor.execute("SELECT  COUNT(*) FROM waiting").fetchone()[0] == 1
    assert work_cursor.execute("SELECT  COUNT(*) FROM complete").fetchone()[0] == 0
    assert work_cursor.execute("SELECT  COUNT(*) FROM proc_comp").fetchone()[0] == 1


def test_workerjob_process_finished(hf, monkeypatch):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    monkeypatch.setattr(rdmcl, "cluster2database", lambda *_: True)
    temp_dir = br.TempDir()
    work_db = temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    rdmcl.WORKER_DB = work_db
    rdmcl.WORKER_OUT = temp_dir.path

    seqbuddy = hf.get_data("cteno_panxs")
    alignment = hf.get_data("cteno_panxs_aln")
    sim_scores = hf.get_data("cteno_sim_scores")
    worker = rdmcl.WorkerJob(seqbuddy, "sql_broker")
    worker.heartbeat.id = 1

    work_con = sqlite3.connect(work_db)
    work_cursor = work_con.cursor()
    work_cursor.execute("INSERT INTO proc_comp (hash, master_id) VALUES (?, ?)", (worker.job_id, worker.heartbeat.id))
    work_con.commit()

    # Run through but don't delete files
    seqs_path = os.path.join(temp_dir.path, "%s.seqs" % worker.job_id)
    seqbuddy.write(seqs_path)
    aln_path = os.path.join(temp_dir.path, "%s.aln" % worker.job_id)
    alignment.write(aln_path)
    graph_path = temp_dir.subfile("%s.graph" % worker.job_id)
    sim_scores.to_csv(graph_path, header=None, index=False)

    assert os.path.isfile(seqs_path)
    assert os.path.isfile(aln_path)
    assert os.path.isfile(graph_path)

    ret_sim_scores, ret_alignment = worker.process_finished()
    assert str(ret_sim_scores) == str(sim_scores)
    assert str(ret_alignment) == str(alignment)

    assert not os.path.isfile(seqs_path)
    assert not os.path.isfile(aln_path)
    assert not os.path.isfile(graph_path)


def test_workerjob_check_if_active(hf, monkeypatch, capsys):
    monkeypatch.setattr(rdmcl, "HeartBeat", MockHeartBeat)
    temp_dir = br.TempDir()

    seqbuddy = hf.get_data("cteno_panxs")
    worker = rdmcl.WorkerJob(seqbuddy, "sql_broker")
    worker.heartbeat.id = 1

    work_db = temp_dir.copy_to("%swork_db.sqlite" % hf.resource_path)
    work_con = sqlite3.connect(work_db)
    work_cursor = work_con.cursor()
    rdmcl.WORKER_DB = work_db
    rdmcl.WORKER_OUT = temp_dir.path

    hb_db = temp_dir.copy_to("%sheartbeat_db.sqlite" % hf.resource_path)
    hb_con = sqlite3.connect(os.path.join(temp_dir.path, "heartbeat_db.sqlite"))
    hb_cursor = hb_con.cursor()
    rdmcl.HEARTBEAT_DB = hb_db

    # No workers no seq file, die
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", (worker.job_id, 1,))
    work_cursor.execute("INSERT INTO waiting (hash, master_id) VALUES (?, ?)", (worker.job_id, 2,))
    work_con.commit()
    assert not worker.check_if_active()
    assert not work_cursor.execute("SELECT * FROM waiting WHERE master_id=1").fetchone()
    assert work_cursor.execute("SELECT * FROM waiting WHERE master_id=2").fetchone()

    # No workers, delete seqs file
    seq_file = temp_dir.subfile("%s.seqs" % worker.job_id)
    worker.heartbeat.id = 2
    assert not worker.check_if_active()
    assert not work_cursor.execute("SELECT * FROM waiting").fetchone()
    assert not os.path.isfile(seq_file)

    # Workers alive, but job is missing
    monkeypatch.setattr(rdmcl.WorkerJob, "pull_from_db", lambda *_: print("Missing job (from monkeypatch)"))
    monkeypatch.setattr(rdmcl.WorkerJob, "queue_job", lambda *_: print("restarting job (from monkeypatch)"))
    hb_cursor.execute("INSERT INTO heartbeat (thread_id, thread_type, pulse) VALUES (?, ?, ?)", (1, "worker", 10e12,))
    hb_con.commit()
    assert worker.check_if_active()
    out, err = capsys.readouterr()
    assert "a2aaca4f79bd56fbf8debfdc281660fd vanished!\nrestarting job (from monkeypatch)\n" in out
    assert "Missing job (from monkeypatch)" in out

    # Job already queued
    work_cursor.execute("INSERT INTO queue (hash) VALUES ('a2aaca4f79bd56fbf8debfdc281660fd')")
    work_con.commit()
    assert worker.check_if_active()
    out, err = capsys.readouterr()
    assert not out

    # Job already complete
    work_cursor.execute("INSERT INTO complete (hash) VALUES ('a2aaca4f79bd56fbf8debfdc281660fd')")
    work_cursor.execute("DELETE FROM queue WHERE hash='a2aaca4f79bd56fbf8debfdc281660fd'")
    work_con.commit()
    assert worker.check_if_active()
    out, err = capsys.readouterr()
    assert not out

    # Confirm a processing job has an active worker on it
    work_cursor.execute("INSERT INTO processing (hash) VALUES ('a2aaca4f79bd56fbf8debfdc281660fd')")
    work_con.commit()
    capsys.readouterr()
    assert worker.check_if_active()
    out, err = capsys.readouterr()
    assert out == "restarting job (from monkeypatch)\n", print(out)

    # Confirm a proc_comp job has an active master on it
    work_cursor.execute("DELETE FROM processing WHERE hash='a2aaca4f79bd56fbf8debfdc281660fd'")
    work_cursor.execute("DELETE FROM complete WHERE hash='a2aaca4f79bd56fbf8debfdc281660fd'")
    work_cursor.execute("INSERT INTO proc_comp (hash, master_id) VALUES ('a2aaca4f79bd56fbf8debfdc281660fd', 5)")
    work_con.commit()
    capsys.readouterr()
    assert worker.check_if_active()
    assert work_cursor.execute("SELECT * FROM complete WHERE hash='a2aaca4f79bd56fbf8debfdc281660fd'").fetchone()
    assert not work_cursor.execute("SELECT * FROM proc_comp WHERE hash='a2aaca4f79bd56fbf8debfdc281660fd'").fetchone()


# #########  MCL stuff  ########## #
def test_mcmcmc_mcl(hf):
    # Need to monkeypatch Cluster, mc_create_all_by_all_scores, Progress.update, helpers.MarkovClustering,
    # and retrieve_all_by_all_scores()

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

    # psi_pred_ss2_paths = hf.get_data("ss2_paths")
    # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqbuddy, psi_pred_ss2_paths, sql_broker)

    cluster = rdmcl.Cluster(cluster_ids, hf.get_db_graph("3c15516819aa19b069b0e8858444f876", sql_broker))
    os.makedirs(os.path.join(ext_tmp_dir.path, "progress"))
    progress = rdmcl.Progress(os.path.join(ext_tmp_dir.path, "progress"), cluster)

    args = (6.372011782427792, 0.901221218627, 1)  # inflation, gq, r_seed
    params = [ext_tmp_dir.path, seqbuddy, cluster, taxa_sep, sql_broker, hf.get_data("ss2_paths"), progress, 3]

    assert rdmcl.mcmcmc_mcl(args, params) == 19.538461538461537
    with open(os.path.join(ext_tmp_dir.path, "max.txt"), "r") as ifile:
        output = ifile.read()
        assert output == "6b39ebc4f5fe7dfef786d8ee3e1594ed,cb23cf3b4d355140e525a1158af5102d," \
                         "8521872b6c07205f3198bb70699f3d93,a06adee8cc3631773890bb5842bf8df9," \
                         "09cac4f034df8a2805171e1e61cc8666", print(output)

    args = (3.1232, 0.73432, 1)  # inflation, gq, r_seed
    assert rdmcl.mcmcmc_mcl(args, params) == 20.923076923076923
    with open(os.path.join(ext_tmp_dir.path, "max.txt"), "r") as ifile:
        output = ifile.read()
        assert output == "6b39ebc4f5fe7dfef786d8ee3e1594ed,cb23cf3b4d355140e525a1158af5102d," \
                         "8521872b6c07205f3198bb70699f3d93,a06adee8cc3631773890bb5842bf8df9," \
                         "09cac4f034df8a2805171e1e61cc8666\n" \
                         "fdffc226205a71dfbdd6cc093cab56cf,cb23cf3b4d355140e525a1158af5102d," \
                         "a06adee8cc3631773890bb5842bf8df9,09cac4f034df8a2805171e1e61cc8666", print(output)

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
    assert hf.string2hash(tmp_file.read()) == "4024a3aa15d605fd8924d8fc6cb9361e", print(tmp_file.read())


# #########  Final sequence placement  ########## #
def test_instantiate_seqs2clusters(hf, monkeypatch):
    # Get starting graph from pre-computed database
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    clusters = hf.get_test_clusters(broker, parent_sb, rdmcl)
    cluster1, cluster2, cluster3, cluster4, cluster5 = clusters

    tmpdir = br.TempDir()
    seq2clust_obj = rdmcl.Seqs2Clusters(clusters, 3, parent_sb, tmpdir.path)

    assert seq2clust_obj.clusters == clusters
    assert seq2clust_obj.min_clust_size == 3
    assert seq2clust_obj.seqbuddy == parent_sb
    assert seq2clust_obj.outdir == tmpdir.path
    assert type(seq2clust_obj.tmp_file) == br.TempFile
    assert type(seq2clust_obj.tmp_dir) == br.TempDir
    assert seq2clust_obj.tmp_dir.subfiles == ["seqs.fa"]
    assert type(seq2clust_obj.small_clusters) == OrderedDict
    assert list(seq2clust_obj.small_clusters.keys()) == ['group_0_2', 'group_0_3', 'group_0_4']
    assert type(seq2clust_obj.large_clusters) == OrderedDict
    assert list(seq2clust_obj.large_clusters.keys()) == ['group_0_0', 'group_0_1']

    small_clusters = [clust.seq_ids for clust_name, clust in seq2clust_obj.small_clusters.items()]
    assert small_clusters == [cluster3.seq_ids, cluster4.seq_ids, cluster5.seq_ids]
    large_clusters = [clust.seq_ids for clust_name, clust in seq2clust_obj.large_clusters.items()]
    assert large_clusters == [cluster1.seq_ids, cluster2.seq_ids]
    assert len([seq_id for sub_clust in small_clusters + large_clusters for seq_id in sub_clust]) == 15
    broker.close()


def test_mc_build_cluster_nulls(hf):
    rsquare_vals_df = pd.read_csv(os.path.join(hf.resource_path, "hmms", "full_r2.csv"), index_col=0)
    global_null_file = br.TempFile()
    cluster_nulls_file = br.TempFile()
    out_of_cluster_file = br.TempFile()
    temp_log_output = br.TempFile()
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    clusters = hf.get_test_clusters(broker, parent_sb, rdmcl)
    broker.close()
    clusters[0].seq_ids = sorted(clusters[0].seq_ids)[:3]
    rsquare_vals_df = rsquare_vals_df.loc[(-rsquare_vals_df["rec_id1"].isin(clusters[1].seq_ids)) &
                                          (-rsquare_vals_df["rec_id2"].isin(clusters[1].seq_ids))]
    parent_sb.records = [rec for rec in parent_sb.records if rec.id not in clusters[1].seq_ids]

    del clusters[1]
    seq2clust_obj = rdmcl.Seqs2Clusters(clusters, 3, parent_sb, br.TempDir().path)

    args = [rsquare_vals_df, global_null_file.path, cluster_nulls_file.path,
            out_of_cluster_file.path, temp_log_output.path]
    seq2clust_obj._mc_build_cluster_nulls(clusters[0], args)

    assert global_null_file.read() == '''\
BOL-PanxαB,Bab-PanxαA,0.9854500607340071
BOL-PanxαB,Bch-PanxαA,0.9814339761739038
BOL-PanxαB,BOL-PanxαB,1.0
Bab-PanxαA,Bch-PanxαA,0.9963598151355412
Bab-PanxαA,Bab-PanxαA,1.0
Bch-PanxαA,Bch-PanxαA,1.0
''', print(global_null_file.read())
    assert cluster_nulls_file.read() == '"group_0_0":{"mu":0.987747950681,"sigma":0.006306366668},', \
        print(cluster_nulls_file.read())
    assert out_of_cluster_file.read() == """\
BOL-PanxαB,Bfo-PanxαE,0.9777341546918392
BOL-PanxαB,Bfr-PanxαA,0.9237171724424902
BOL-PanxαB,Hca-PanxαA,0.97976788883183
BOL-PanxαB,Hvu-PanxβA,0.2369553465136797
BOL-PanxαB,Lcr-PanxαG,0.9864552763789216
BOL-PanxαB,Vpa-PanxαD,0.07725859240096591
BOL-PanxαB,Oma-PanxαB,0.03632897706786227
Bab-PanxαA,Bfo-PanxαE,0.977670073270024
Bab-PanxαA,Bfr-PanxαA,0.9454631199836152
Bab-PanxαA,Hca-PanxαA,0.994053564420282
Bab-PanxαA,Hvu-PanxβA,0.2258542797590169
Bab-PanxαA,Lcr-PanxαG,0.9990117173041204
Bab-PanxαA,Vpa-PanxαD,0.0852327674206815
Bab-PanxαA,Oma-PanxαB,0.039921601139418186
Bch-PanxαA,Bfo-PanxαE,0.9750441560758988
Bch-PanxαA,Bfr-PanxαA,0.9367784583942784
Bch-PanxαA,Hca-PanxαA,0.9921229245873938
Bch-PanxαA,Hvu-PanxβA,0.22872882053802146
Bch-PanxαA,Lcr-PanxαG,0.9969676816010246
Bch-PanxαA,Vpa-PanxαD,0.08142126155716654
Bch-PanxαA,Oma-PanxαB,0.03978017203455591
""", print(out_of_cluster_file.read())
    assert temp_log_output.read() == """\
group_0_0
\tN: 3
\tMean: 0.987747950681
\tStd: 0.006306366668

""", print(temp_log_output.read())

    temp_log_output.clear()
    seq2clust_obj._mc_build_cluster_nulls(clusters[1], args)
    assert temp_log_output.read() == """\
group_0_2
\tN: 1
\tMean: Null
\tStd: Null

""", print(temp_log_output.read())


def test_mc_build_seq2group(hf):
    rsquare_vals_df = pd.read_csv(os.path.join(hf.resource_path, "hmms", "full_r2.csv"), index_col=0)
    seq2group_dists_file = br.TempFile()
    orig_clusters_file = br.TempFile()
    temp_log_output = br.TempFile()
    args = [rsquare_vals_df, seq2group_dists_file.path, orig_clusters_file.path, temp_log_output.path]
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    clusters = hf.get_test_clusters(broker, parent_sb, rdmcl)
    broker.close()
    seq2clust_obj = rdmcl.Seqs2Clusters(clusters, 3, parent_sb, br.TempDir().path)

    seq2clust_obj._mc_build_seq2group("Bab-PanxαA", args)

    assert seq2group_dists_file.read() == '''\
"Bab-PanxαA":{"group_0_0":0.788972934053,"group_0_1":0.074113151245,"group_0_2":0.085232767421,\
"group_0_3":0.996532640862,"group_0_4":0.225854279759},''', print(seq2group_dists_file.read())
    assert orig_clusters_file.read() == '"Bab-PanxαA":"group_0_0",'
    assert temp_log_output.read() == """\
Bab-PanxαA
\tgroup_0_0: 0.788972934053
\tgroup_0_1: 0.074113151245
\tgroup_0_2: 0.085232767421
\tgroup_0_3: 0.996532640862
\tgroup_0_4: 0.225854279759
\tBest: group_0_3
""", print(temp_log_output.read())


def test_create_hmms_for_every_rec(hf):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    clusters = hf.get_test_clusters(broker, parent_sb, rdmcl)

    tmpdir = br.TempDir()
    hmm_dir = tmpdir.subdir("hmm")
    with open(os.path.join(hmm_dir, "Bab-PanxαA.hmm"), "w") as ofile:
        ofile.write("Testing testing, 1, 2, 3...")
    seq2clust_obj = rdmcl.Seqs2Clusters(clusters, 3, parent_sb, tmpdir.path)
    seq2clust_obj.seqbuddy.records = seq2clust_obj.seqbuddy.records[0:2]
    seq2clust_obj.create_hmms_for_every_rec()
    assert sorted(os.listdir(os.path.join(seq2clust_obj.outdir, "hmm"))) == ['BOL-PanxαB.hmm', 'Bab-PanxαA.hmm']
    with open(os.path.join(hmm_dir, "Bab-PanxαA.hmm"), "r") as ifile:
        assert ifile.read() == "Testing testing, 1, 2, 3..."
    with open(os.path.join(hmm_dir, "BOL-PanxαB.hmm"), "r") as ifile:
        assert "HMM          A        C        D        E        F        G        H        I        K" in ifile.read()
    broker.close()


def test_create_hmm_fwd_score_df(hf, monkeypatch):
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    parent_sb.records = [rec for rec in parent_sb.records if rec.id in
                         ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE']]
    seq2clust_obj = rdmcl.Seqs2Clusters([], 3, parent_sb, "foo.path")

    monkeypatch.setattr(rdmcl.Seqs2Clusters, "_separate_large_small", lambda *_: True)
    monkeypatch.setattr(rdmcl.Seqs2Clusters, "create_hmms_for_every_rec",
                        lambda *_: os.path.join(hf.resource_path, "hmms"))
    hmm_fwd_scores = seq2clust_obj.create_hmm_fwd_score_df()
    # hmm_fwd_scores.to_csv("tests/unit_test_resources/hmms/fwd_df.csv")
    assert str(hmm_fwd_scores) == """\
        hmm_id      rec_id   fwd_raw
0   BOL-PanxαB  BOL-PanxαB  702.0341
1   BOL-PanxαB  Bab-PanxαA  641.0910
2   BOL-PanxαB  Bch-PanxαA  626.8386
3   BOL-PanxαB  Bfo-PanxαE  619.8560
4   Bab-PanxαA  BOL-PanxαB  636.7927
5   Bab-PanxαA  Bab-PanxαA  696.4839
6   Bab-PanxαA  Bch-PanxαA  663.5856
7   Bab-PanxαA  Bfo-PanxαE  623.7233
8   Bch-PanxαA  BOL-PanxαB  627.4737
9   Bch-PanxαA  Bab-PanxαA  668.5334
10  Bch-PanxαA  Bch-PanxαA  702.3412
11  Bch-PanxαA  Bfo-PanxαE  619.3785
12  Bfo-PanxαE  BOL-PanxαB  615.9315
13  Bfo-PanxαE  Bab-PanxαA  624.0910
14  Bfo-PanxαE  Bch-PanxαA  615.3356
15  Bfo-PanxαE  Bfo-PanxαE  702.2801""", print(str(hmm_fwd_scores))


def test_create_fwd_score_rsquared_matrix(hf, monkeypatch):
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    parent_sb.records = [rec for rec in parent_sb.records if rec.id in
                         ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE']]
    seq2clust_obj = rdmcl.Seqs2Clusters([], 3, parent_sb, "foo.path")
    hmm_fwd_scores = pd.read_csv(os.path.join(hf.resource_path, "hmms", "fwd_df.csv"), index_col=0)

    monkeypatch.setattr(rdmcl.Seqs2Clusters, "_separate_large_small", lambda *_: True)
    monkeypatch.setattr(rdmcl.Seqs2Clusters, "create_hmm_fwd_score_df",
                        lambda *_: hmm_fwd_scores)

    rsquare_vals_df = seq2clust_obj.create_fwd_score_rsquared_matrix()
    # rsquare_vals_df.to_csv("tests/unit_test_resources/hmms/fwd_r2.csv")
    assert str(rsquare_vals_df) == """\
      rec_id1     rec_id2        r_square
0  BOL-PanxαB  Bab-PanxαA  0.016894041431
1  BOL-PanxαB  Bch-PanxαA  0.087311057754
2  BOL-PanxαB  Bfo-PanxαE  0.274041115357
3  BOL-PanxαB  BOL-PanxαB  1.000000000000
4  Bab-PanxαA  Bch-PanxαA  0.497451379268
5  Bab-PanxαA  Bfo-PanxαE  0.453877689738
6  Bab-PanxαA  Bab-PanxαA  1.000000000000
7  Bch-PanxαA  Bfo-PanxαE  0.390838714109
8  Bch-PanxαA  Bch-PanxαA  1.000000000000
9  Bfo-PanxαE  Bfo-PanxαE  1.000000000000""", print(rsquare_vals_df)


def test_create_truncnorm(hf):
    hmm_fwd_scores = pd.read_csv(os.path.join(hf.resource_path, "hmms", "fwd_r2.csv"), index_col=0)
    truncnorm = rdmcl.Seqs2Clusters._create_truncnorm(helpers.mean(hmm_fwd_scores.r_square),
                                                      helpers.std(hmm_fwd_scores.r_square))
    assert truncnorm.pdf(0.8) == 1.0920223528477309


def test_place_seqs_in_clusts(hf, monkeypatch):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)
    clusters = hf.get_test_clusters(broker, parent_sb, rdmcl)
    subset_ids = [seq for clust in clusters for seq in clust.seq_ids]
    parent_sb.records = [rec for rec in parent_sb.records if rec.id in subset_ids]
    tmpdir = br.TempDir()
    tmpdir.subdir("hmm")
    seq2clust_obj = rdmcl.Seqs2Clusters(clusters, 3, parent_sb, tmpdir.path)

    hmm_fwd_scores = pd.read_csv(os.path.join(hf.resource_path, "hmms", "full_r2.csv"), index_col=0)
    monkeypatch.setattr(rdmcl.Seqs2Clusters, "create_fwd_score_rsquared_matrix",
                        lambda *_: hmm_fwd_scores)

    seq2clust_obj.place_seqs_in_clusts()
    all_clusts = sorted([sorted(c.seq_ids) for c in seq2clust_obj.clusters])
    assert len(seq2clust_obj.clusters) == 4, print(all_clusts)
    cstring = "\n".join([str(c) for c in all_clusts])
    assert all_clusts[0] == ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE',
                             'Bfr-PanxαA', 'Hca-PanxαA', 'Lcr-PanxαG'], print(cstring)
    assert all_clusts[1] == ["Hvu-PanxβA"], print(cstring)
    assert all_clusts[2] == ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD',
                             'Pba-PanxαB', 'Tin-PanxαF', 'Vpa-PanxαD'], print(cstring)
    assert all_clusts[3] == ["Oma-PanxαB"], print(cstring)

    # assert hf.string2hash(seq2clust_obj.tmp_file.read()) == "92d4b9600d74732452c1d81b7d1a8ece"


def test_get_sort_order():
    seq2group_dists = {"Bab-PanxαC": {"group_0_0": 0.088972934053, "group_0_1": 0.074113151245,
                                      "group_0_2": 0.085232767421, "group_0_3": 0.496532640862,
                                      "group_0_4": 0.935854279759},
                       "Bab-PanxαA": {"group_0_0": 0.788972934053, "group_0_1": 0.074113151245,
                                      "group_0_2": 0.085232767421, "group_0_3": 0.996532640862,
                                      "group_0_4": 0.225854279759},
                       "Bab-PanxαB": {"group_0_0": 0.088972934053, "group_0_1": 0.074113151245,
                                      "group_0_2": 0.085232767421, "group_0_3": 0.496532640862,
                                      "group_0_4": 0.925854279759}}
    new_group = {"group_0_0": {"Mle-PanxαA", "Tin-PanxαA", "Dgl-PanxαA"},
                 "group_0_1": {"Mle-PanxαC", "Tin-PanxαC", "Dgl-PanxαC"},
                 "group_0_2": {"Mle-PanxαD", "Tin-PanxαD", "Dgl-PanxαD"},
                 "group_0_3": {"Bab-PanxαA", "Bab-PanxαC", "Tin-PanxαE", "Dgl-PanxαE"},
                 "group_0_4": {"Bab-PanxαB", "Tin-PanxαE", "Dgl-PanxαE"}}
    placed = []
    orig_clusts = {"Bab-PanxαA": "group_0_3", "Bab-PanxαB": "group_0_4", "Bab-PanxαC": "group_0_3"}
    sort_order = rdmcl.Seqs2Clusters.get_sort_order(seq2group_dists, placed, new_group, orig_clusts)
    assert sort_order == [[('Bab-PanxαB', 0.925854279759)], [('Bab-PanxαA', 0.996532640862),
                                                             ('Bab-PanxαC', 0.935854279759)]], print(sort_order)

    placed = ['Bab-PanxαA']
    sort_order = rdmcl.Seqs2Clusters.get_sort_order(seq2group_dists, placed, new_group, orig_clusts)
    assert sort_order == [[('Bab-PanxαB', 0.925854279759)], [('Bab-PanxαC', 0.935854279759)]], print(sort_order)


def test_place_orphans(hf):
    return
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()

    graph = hf.get_db_graph("6935966a6b9967c3006785488d968230", broker)
    parent_sb = rdmcl.Sb.SeqBuddy("%sCteno_pannexins_subset.fa" % hf.resource_path)

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

    orphans = rdmcl.Orphans(parent_sb, clusters, broker, hf.get_data("ss2_paths"), br.TempDir().path)

    # ##### mc_check_orphans() ##### #
    tmp_file = br.TempFile()
    orphans.mc_check_orphans(cluster3, [tmp_file.path])
    orphans.mc_check_orphans(cluster4, [tmp_file.path])
    orphans.mc_check_orphans(cluster5, [tmp_file.path])
    assert tmp_file.read() == """\
group_0_2\tgroup_0_1\t0.0331066411881
group_0_3\tgroup_0_0\t0.0174886528987
""", print(tmp_file.read())

    assert hf.string2hash(orphans.tmp_file.read()) == "b57ba819b98f110758c4abc31c36a906", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### _check_orphan() ##### #
    # Cluster is placed
    orphs = orphans._check_orphan(cluster3)

    assert orphs[0] == 'group_0_1'
    assert round(orphs[1], 12) == 0.033106641188
    assert not orphans._check_orphan(cluster5)  # Insufficient support for largest cluster
    assert hf.string2hash(orphans.tmp_file.read()) == "625759dd2017515cad4c9d1200311fd2", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### place_orphans() ##### #
    orphans.place_orphans()
    assert orphans.clusters[0].seq_ids == ["Hvu-PanxβA"]
    assert orphans.clusters[1].seq_ids == ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE',
                                           'Bfr-PanxαA', 'Hca-PanxαA', 'Lcr-PanxαG']
    assert orphans.clusters[2].seq_ids == ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD',
                                           'Pba-PanxαB', 'Tin-PanxαF', 'Vpa-PanxαD']
    assert hf.string2hash(orphans.tmp_file.read()) == "92d4b9600d74732452c1d81b7d1a8ece", print(orphans.tmp_file.read())

    # Multicore doesn't seem to work in py.test, but can at least call it like it does
    clusters = [deepcopy(clust) for clust in orig_clusters]
    orphans = rdmcl.Orphans(parent_sb, clusters, broker, hf.get_data("ss2_paths"), br.TempDir().path)
    orphans.place_orphans()
    assert len(orphans.clusters) == 3

    # If no small clusters, nothing much happens
    orphans = rdmcl.Orphans(parent_sb, clusters[:2], broker, hf.get_data("ss2_paths"), br.TempDir().path)
    assert not orphans.place_orphans()
    assert len(orphans.clusters) == 2


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
parser.add_argument("-r", "--resume", action="store_true",
                    help="Try to pick up where a previous run left off (this breaks r_seed).")
parser.add_argument("-trm", "--trimal", action="append", nargs="+",
                    help="Specify a list of trimal thresholds to apply (move from more strict to less)")
parser.add_argument("-f", "--force", action="store_true",
                    help="Overwrite previous run")
parser.add_argument("-q", "--quiet", action="store_true",
                    help="Suppress all output during run (only final output is returned)")

# This is to allow py.test to work with its own flags
in_args = parser.parse_args([])


def test_argparse_init(monkeypatch, hf):
    out_dir = br.TempDir()
    argv = ['rdmcl.py', os.path.join(hf.resource_path, "BOL_Lcr_Mle_Vpa.fa"), "-o", out_dir.path]
    monkeypatch.setattr(rdmcl.sys, "argv", argv)
    temp_in_args = rdmcl.argparse_init()
    assert temp_in_args.mcmc_steps == 0
    assert temp_in_args.open_penalty == -5
    assert temp_in_args.ext_penalty == 0


@pytest.mark.slow
def test_full_run(hf, capsys):
    # I can't break these up into separate test functions because of collisions with logger
    out_dir = br.TempDir()

    # First try a RESUME run
    test_in_args = deepcopy(in_args)
    test_in_args.sequences = os.path.join(hf.resource_path, "BOL_Lcr_Mle_Vpa.fa")
    test_in_args.outdir = out_dir.path
    test_in_args.sqlite_db = os.path.join(hf.resource_path, "db.sqlite")
    test_in_args.psipred_dir = os.path.join(hf.resource_path, "psi_pred")
    test_in_args.mcmc_steps = 10
    test_in_args.r_seed = 1
    rdmcl.full_run(test_in_args)

    for expected_dir in ["alignments", "hmm", "mcmcmc", "psi_pred", "sim_scores"]:
        expected_dir = os.path.join(out_dir.path, expected_dir)
        assert os.path.isdir(expected_dir), print(expected_dir)

    for expected_file in ["final_clusters.txt", "placement.log", "paralog_cliques", "rdmcl.log"]:
        expected_file = os.path.join(out_dir.path, expected_file)
        assert os.path.isfile(expected_file), print(expected_file)

    with open(os.path.join(out_dir.path, "final_clusters.txt"), "r") as ifile:
        content = ifile.read()
        assert content == """\
group_0_0_0\t20.3333\tBOL-PanxαA\tLcr-PanxαH\tMle-Panxα10A\tMle-Panxα9\tVpa-PanxαB
group_0_2\t4.0\tBOL-PanxαH\tMle-Panxα8
group_0_1_1\t5.25\tLcr-PanxαK\tMle-Panxα7A
group_0_0_1\t1.25\tMle-Panxα5
group_0_1_0\t2.0833\tBOL-PanxαG
""", print(content)

    out, err = capsys.readouterr()
    print(err)
    assert "RESUME: All PSI-Pred .ss2 files found" in err
    assert "Generating initial all-by-all similarity graph" in err

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
    subset_ids = ['BOL-PanxαA', 'BOL-PanxαG', 'BOL-PanxαH', 'Lcr-PanxαH', 'Lcr-PanxαK', 'Mle-Panxα5',
                  'Mle-Panxα7A', 'Mle-Panxα8', 'Mle-Panxα9', 'Mle-Panxα10A', 'Vpa-PanxαB']
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

    # shutil.copyfile(os.path.join(out_dir.path, "sqlite_db.sqlite"), "unit_test_db")

    for expected_dir in ['alignments', 'hmm', 'mcmcmc', 'psi_pred', 'sim_scores']:
        assert os.path.isdir(os.path.join(out_dir.path, expected_dir)), print(next(os.walk(out_dir.path))[1])

    for expected_file in ["final_clusters.txt", "orphans.log",
                          "paralog_cliques", "rdmcl.log", "sqlite_db.sqlite"]:
        assert os.path.isfile(os.path.join(out_dir.path, expected_file)), print(next(os.walk(out_dir.path))[2])

    with open(os.path.join(out_dir.path, "final_clusters.txt"), "r") as ifile:
        content = ifile.read()
        assert content == """\
group_0_0_0\t20.3333\tBOL-PanxαA\tLcr-PanxαH\tMle-Panxα10A\tMle-Panxα9\tVpa-PanxαB
group_0_2\t4.0\tBOL-PanxαH\tMle-Panxα8
group_0_1\t9.0417\tLcr-PanxαK\tMle-Panxα7A
group_0_0_1\t1.25\tMle-Panxα5
group_0_3\t2.0833\tBOL-PanxαG
""", print(content)

    out, err = capsys.readouterr()
    assert "Generating initial all-by-all similarity graph (55 comparisons)" in err, print(err)

    # Cleanup
    if os.path.isfile("rdmcl.log"):
        os.remove("rdmcl.log")
    '''
