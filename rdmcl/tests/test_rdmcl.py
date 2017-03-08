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
from copy import copy
import shutil


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
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    assert [taxa for taxa in cluster.taxa] == ['BOL', 'Bab', 'Bch', 'Bfo', 'Bfr', 'Cfu', 'Dgl', 'Edu', 'Hca', 'Hru',
                                               'Hvu', 'Lcr', 'Lla', 'Mle', 'Oma', 'Pba', 'Tin', 'Vpa']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "984c4424c2b8529694696d715c4108a5"
    assert cluster.taxa_separator == "-"
    assert cluster.parent is None
    assert cluster.subgroup_counter == 0
    assert cluster.cluster_score is None
    assert cluster.collapsed_genes == OrderedDict([('Cfu-PanxαA', ['Cfu-PanxαF']),
                                                   ('Hvu-PanxβM', ['Hvu-PanxβI', 'Hvu-PanxβG', 'Hvu-PanxβH',
                                                                   'Hvu-PanxβD', 'Hvu-PanxβK', 'Hvu-PanxβE',
                                                                   'Hvu-PanxβF', 'Hvu-PanxβA', 'Hvu-PanxβC',
                                                                   'Hvu-PanxβB', 'Hvu-PanxβJ', 'Hvu-PanxβL',
                                                                   'Hvu-PanxβO']), ('Lcr-PanxαA', ['Lcr-PanxαL']),
                                                   ('Mle-Panxα10A', ['Mle-Panxα9'])])
    assert cluster._name == "group_0"
    assert cluster.seq_ids == ['BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαE', 'BOL-PanxαF',
                               'BOL-PanxαG', 'BOL-PanxαH', 'Bab-PanxαA', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαD',
                               'Bab-PanxαE', 'Bch-PanxαA', 'Bch-PanxαB', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαE',
                               'Bfo-PanxαA', 'Bfo-PanxαB', 'Bfo-PanxαC', 'Bfo-PanxαD', 'Bfo-PanxαE', 'Bfo-PanxαF',
                               'Bfo-PanxαG', 'Bfo-PanxαH', 'Bfo-PanxαI', 'Bfo-PanxαJ', 'Bfr-PanxαA', 'Bfr-PanxαB',
                               'Bfr-PanxαC', 'Bfr-PanxαD', 'Cfu-PanxαA', 'Cfu-PanxαB', 'Cfu-PanxαC', 'Cfu-PanxαD',
                               'Cfu-PanxαE', 'Dgl-PanxαA', 'Dgl-PanxαB', 'Dgl-PanxαC', 'Dgl-PanxαD', 'Dgl-PanxαE',
                               'Dgl-PanxαF', 'Dgl-PanxαG', 'Dgl-PanxαH', 'Dgl-PanxαI', 'Edu-PanxαA', 'Edu-PanxαB',
                               'Edu-PanxαC', 'Edu-PanxαD', 'Edu-PanxαE', 'Edu-PanxαF', 'Edu-PanxαG', 'Edu-PanxαH',
                               'Edu-PanxαI', 'Hca-PanxαA', 'Hca-PanxαB', 'Hca-PanxαC', 'Hca-PanxαD', 'Hca-PanxαE',
                               'Hca-PanxαF', 'Hca-PanxαG', 'Hca-PanxαH', 'Hru-PanxαA', 'Hru-PanxαB', 'Hru-PanxαC',
                               'Hru-PanxαD', 'Hru-PanxαE', 'Hvu-PanxβM', 'Lcr-PanxαA', 'Lcr-PanxαB', 'Lcr-PanxαC',
                               'Lcr-PanxαD', 'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαG', 'Lcr-PanxαH', 'Lcr-PanxαI',
                               'Lcr-PanxαJ', 'Lcr-PanxαK', 'Lla-PanxαA', 'Lla-PanxαB', 'Lla-PanxαC', 'Mle-Panxα1',
                               'Mle-Panxα10A', 'Mle-Panxα11', 'Mle-Panxα12', 'Mle-Panxα2', 'Mle-Panxα3', 'Mle-Panxα4',
                               'Mle-Panxα5', 'Mle-Panxα6', 'Mle-Panxα7A', 'Mle-Panxα8', 'Oma-PanxαA', 'Oma-PanxαB',
                               'Oma-PanxαC', 'Oma-PanxαD', 'Pba-PanxαA', 'Pba-PanxαB', 'Pba-PanxαC', 'Pba-PanxαD',
                               'Pba-PanxαE', 'Pba-PanxαF', 'Pba-PanxαG', 'Tin-PanxαA', 'Tin-PanxαB', 'Tin-PanxαC',
                               'Tin-PanxαD', 'Tin-PanxαE', 'Tin-PanxαF', 'Vpa-PanxαA', 'Vpa-PanxαB', 'Vpa-PanxαC',
                               'Vpa-PanxαD', 'Vpa-PanxαE', 'Vpa-PanxαF', 'Vpa-PanxαG']
    assert cluster.seq_id_hash == "b7636e7de3a0a96e2631db0ba01c0ffc"

    with pytest.raises(ValueError) as err:
        sim_scores = hf.get_data("cteno_sim_scores")
        sim_scores = sim_scores.ix[1:, :]
        rdmcl.Cluster(cluster.seq_ids, sim_scores)
    assert "The number of incoming sequence ids (118) does not match the expected graph size of 6903" in str(err)


def test_cluster_instantiate_child(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    sim_scores = pd.read_csv("%sCteno_pannexins_subgroup_sim.scores" % hf.resource_path, index_col=False, header=None)
    sim_scores.columns = ["seq1", "seq2", "score"]
    cluster = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert [taxa for taxa in cluster.taxa] == ['BOL', 'Bab', 'Bch', 'Bfo', 'Dgl', 'Edu', 'Hca',
                                               'Hru', 'Lcr', 'Mle', 'Oma', 'Tin', 'Vpa']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "eb1ef296e9f4e74cd6deea490a447326"
    assert cluster.taxa_separator == "-"
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
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    global_best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
    best_hits = cluster.recursive_best_hits('Bab-PanxαB', global_best_hits, ['Bab-PanxαB'])
    assert best_hits.to_csv() == """\
,seq1,seq2,score
0,Bab-PanxαB,Vpa-PanxαB,0.9715263513449742
1,Lcr-PanxαH,Vpa-PanxαB,0.979692672624647
2,Lcr-PanxαH,Vpa-PanxαB,0.979692672624647
"""


def test_cluster_perterb(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    global_best_hits = pd.DataFrame(columns=["seq1", "seq2", "score"])
    best_hits = cluster.recursive_best_hits('Lcr-PanxαH', global_best_hits, ['Lcr-PanxαH'])
    assert best_hits.iloc[0].score == 0.979692672624647

    best_hits = cluster.perturb(best_hits)
    assert best_hits.iloc[0].score != 0.979692672624647
    assert round(best_hits.iloc[0].score, 5) == 0.97969


def test_cluster_get_base_cluster(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())

    # No paralogs
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    # First time calling Cluster.score() calculates the score
    assert child.score() == 61.88340444450972

    # The second call just retrieves the attribute from the cluster saved during first call
    assert child.score() == 61.88340444450972

    # With paralogs
    child_ids = ['BOL-PanxαA', 'BOL-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
                 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child.score() == 44.05560826220814

    # Single sequence
    child_ids = ['BOL-PanxαA']
    child = rdmcl.Cluster(child_ids, hf.get_sim_scores(child_ids), parent=parent)
    assert child.score() == -6.000949946302431

    # Include an orphan sequence
    child_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC']
    sim_scores = hf.get_sim_scores(child_ids)
    child = rdmcl.Cluster(child_ids, sim_scores, parent=parent)
    assert child.score() == 0.5466386920733747
    child.seq_ids.append("Foo-Bar3")
    assert child.score(force=True) == 1.3326109193908637

    # Edge case where child is full size of parent
    child = rdmcl.Cluster(parent.seq_ids, parent.sim_scores, parent=parent)
    assert child.score() == -208.64552691255088


def test_cluster_pull_scores_subgraph(hf):
    clust = rdmcl.Cluster(*hf.base_cluster_args())
    assert str(clust.pull_scores_subgraph(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC'])) == """\
            seq1        seq2     score
132   Hru-PanxαA  Lcr-PanxαH  0.966569
165   Hru-PanxαA  Tin-PanxαC  0.958389
6851  Lcr-PanxαH  Tin-PanxαC  0.974295"""


def test_cluster_create_rbh_cliques(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    log_file = br.TempFile()

    # Cluster of less than 6 seqs returns self
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_0 ####### #
['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC']
\tTERMINATED: Group too small."""
    log_file.clear()

    # If clique brings in entire cluster, return self
    seq_ids = ['Bab-PanxαC', 'Bfo-PanxαC', 'BOL-PanxαF', 'Vpa-PanxαD', 'BOL-PanxαE', 'Mle-Panxα11']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_1 ####### #
['BOL-PanxαE', 'BOL-PanxαF', 'Bab-PanxαC', 'Bfo-PanxαC', 'Mle-Panxα11', 'Vpa-PanxαD']
\tTERMINATED: Entire cluster pulled into clique on BOL."""
    log_file.clear()

    # If no paralogs in cluster, return self
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE', 'Vpa-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_2 ####### #
['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαE']
\tTERMINATED: No paralogs present."""
    log_file.clear()

    # When all cliques are disqualified due to size and/or overlap
    seq_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB',
               'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα9', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert log_file.read().strip() == """\
# ####### Testing group_0_3 ####### #
['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH', 'Mle-Panxα10A', 'Mle-Panxα9', 'Oma-PanxαC', 'Tin-PanxαC', 'Vpa-PanxαB']
\tChecking taxa for overlapping cliques:
\t\t# #### Mle #### #
\t\tDisqualified cliques:
\t\t\t['Mle-Panxα9', 'Vpa-PanxαB'] is too small
\t\t\t['Mle-Panxα10A', 'Mle-Panxα9', 'Vpa-PanxαB'] overlaps with ['Mle-Panxα9', 'Vpa-PanxαB']

\t\t!! ALL CLIQUES DISQUALIFIED !!

\tTERMINATED: No significant cliques identified."""
    log_file.clear()

    # When spinning off an otherwise valid clique would orphan one or more sequences, return self
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI', 'Bch-PanxαB',  'Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    assert cluster.create_rbh_cliques(log_file)[0] == cluster
    assert "TERMINATED: Spinning off cliques would orphan Bch-PanxαB." in log_file.read()

    # Check whether identified cliques can be combined
    seq_ids = ['Vpa-PanxαA', 'Mle-Panxα4', 'BOL-PanxαF', 'Lcr-PanxαI', 'Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids), parent=parent)
    cluster.set_name()
    cliques = cluster.create_rbh_cliques(log_file)
    assert cliques[0].seq_ids == ['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4', 'Vpa-PanxαA']
    assert cliques[1].seq_ids == ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']
    assert """Checking if new cliques can combine
\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4'] - ['Lcr-PanxαL', 'Mle-Panxα5', 'Vpa-PanxαF']  NO
\t\t['BOL-PanxαF', 'Lcr-PanxαI', 'Mle-Panxα4'] - ['BOL-PanxαF', 'Mle-Panxα4', 'Vpa-PanxαA']  YES""" in log_file.read()


def test_cluster_len(hf):
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids))
    assert len(cluster) == 5


def test_cluster_str(hf):
    seq_ids = ['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']
    cluster = rdmcl.Cluster(seq_ids, hf.get_sim_scores(seq_ids))
    assert str(cluster) == "['Dgl-PanxαE', 'Hru-PanxαA', 'Lcr-PanxαH', 'Oma-PanxαC', 'Tin-PanxαC']"


def test_cluster2database(hf):
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
    assert response[0][0] == "b7636e7de3a0a96e2631db0ba01c0ffc"                 # hash
    assert 'BOL-PanxαA, BOL-PanxαB, BOL-PanxαC, BOL-PanxαD' in response[0][1]   # seq_ids
    assert response[0][2] == '>Seq1\nMPQQCS-SS\n>Seq2\nMPQICMAAS'               # alignment
    assert 'Hca-PanxαG,Lla-PanxαC,0.42864589074736\n' in response[0][3]         # graph
    # assert response[0][4] == '-55324684799997.984'  # ToDo: decide on final scoring system, or remove  # score
    connect.close()

# #########  PSI-PRED  ########## #
bins = ['chkparse', 'psipass2', 'psipred', 'seq2mtx']


@pytest.mark.parametrize("binary", bins)
def test_psipred_bins(binary, hf):
    assert os.path.isfile("{0}..{1}..{1}psipred{1}bin{1}{2}".format(hf.resource_path, hf.sep, binary))


def test_run_psi_pred(hf):
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Bab-PanxαA")
    ss2_file = rdmcl.run_psi_pred(seqbuddy.records[0])
    assert hf.string2hash(ss2_file) == "b50f39dc22e4d16be325efdd14f7900d"


def test_mc_psi_pred(hf):
    tmpdir = br.TempDir()
    tmpdir.subdir("psi_pred")
    rdmcl.mc_psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}psi_pred{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.mc_psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}psi_pred{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
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


# #########  Miscellaneous  ########## #
def test_progress(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    tmpdir = br.TempDir()
    progress = rdmcl.Progress(tmpdir.path, cluster)
    assert progress.outdir == tmpdir.path
    assert os.path.isfile("{0}{1}.progress".format(tmpdir.path, hf.sep))
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # The dictionary is not static, so just sort the string:
        # {"placed": 0, "mcl_runs": 0, "total": 118}
        assert "".join(sorted(ifile.read())) == '     """""",,00118:::_aaccdelllmnoprsttu{}'

    progress.update("mcl_runs", 2)
    with open("{0}{1}.progress".format(tmpdir.path, hf.sep), "r") as ifile:
        # {"placed": 0, "mcl_runs": 2, "total": 118}
        assert "".join(sorted(ifile.read())) == '     """""",,01128:::_aaccdelllmnoprsttu{}'

    json = progress.read()
    assert json["mcl_runs"] == 2
    assert json["placed"] == 0
    assert json["total"] == 118

    assert str(progress) == "MCL runs processed: 2. Sequences placed: 0/118. Run time: "


def test_check_sequences(hf, monkeypatch, capsys):
    monkeypatch.setattr(rdmcl, "logging", MockLogging)
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.check_sequences(seqbuddy, "-")
    out, err = capsys.readouterr()
    assert "Checking that the format of all sequence ids matches 'taxa-gene'" in out
    assert "    134 sequences PASSED" in out

    rdmcl.Sb.rename(seqbuddy, "Mle-Panxα1", "Mle:Panxα1")
    with pytest.raises(SystemExit):
        rdmcl.check_sequences(seqbuddy, "-")

    out, err = capsys.readouterr()
    assert "Checking that the format of all sequence ids matches 'taxa-gene'" in out
    print(out)
    assert "Malformed sequence id(s): 'Mle:Panxα1, Mle:Panxα10A, Mle:Panxα11, Mle:Panxα12'\n" \
           "The taxa separator character is currently set to '-',\n" \
           " which can be changed with the '-ts' flag" in out


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


def test_mc_score_sequence(hf):
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    alignbuddy = rdmcl.Alb.generate_msa(seqbuddy, "mafft", params="--globalpair --thread -2", quiet=True)
    psi_pred_files = [(rec.id, rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, hf.sep, rec.id)))
                      for rec in alignbuddy.records()]
    psi_pred_files = OrderedDict(psi_pred_files)
    tmp_dir = br.TempDir()

    seq_pairs = [("Mle-Panxα2", "Mle-Panxα12"), ("Mle-Panxα1", "Mle-Panxα3")]
    args = [alignbuddy, psi_pred_files, tmp_dir.path, -5, 0]
    rdmcl.mc_score_sequences(seq_pairs, args)
    assert os.path.isfile("%s%s50687872cdaaed1fee7e809b5a032377" % (tmp_dir.path, hf.sep))
    with open("%s%s50687872cdaaed1fee7e809b5a032377" % (tmp_dir.path, hf.sep), "r") as ifile:
        assert ifile.read() == """
Mle-Panxα2,Mle-Panxα12,0.4789936469710511
Mle-Panxα1,Mle-Panxα3,0.3706607994410156"""


def test_compare_pairwise_alignment():
    seqbuddy = rdmcl.Sb.SeqBuddy(">seq1\nMPQMSASWI\n>Seq2\nMPPQISASI")
    alignbuddy = rdmcl.Alb.generate_msa(seqbuddy, "mafft", params="--globalpair --op 0 --thread -2", quiet=True)
    assert str(alignbuddy) == """\
>seq1
MP-QMSASWI
>Seq2
MPPQISAS-I
"""
    subs_mat_score = rdmcl.compare_pairwise_alignment(alignbuddy, -5, 0)
    assert subs_mat_score == 0.5176252319109462


def test_create_all_by_all_scores(hf):
    seqbuddy = rdmcl.Sb.SeqBuddy(hf.get_data("cteno_panxs"))
    rdmcl.Sb.pull_recs(seqbuddy, "Mle")
    alignbuddy = rdmcl.Alb.generate_msa(seqbuddy, "mafft", params="--globalpair --thread -2", quiet=True)
    psi_pred_files = [(rec.id, rdmcl.read_ss2_file("%spsi_pred%s%s.ss2" % (hf.resource_path, hf.sep, rec.id)))
                      for rec in alignbuddy.records()]
    psi_pred_files = OrderedDict(psi_pred_files)

    sim_scores = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_files)
    assert len(sim_scores.index) == 66  # This is for 12 starting sequences --> (a * (a - 1)) / 2
    compare = sim_scores.loc[:][(sim_scores['seq1'] == "Mle-Panxα2") & (sim_scores['seq2'] == "Mle-Panxα12")]
    assert compare.iloc[0]['score'] == 0.54706454433078899

    rdmcl.Alb.pull_records(alignbuddy, "Mle-Panxα2")
    sim_scores = rdmcl.create_all_by_all_scores(alignbuddy, psi_pred_files)
    assert len(sim_scores.index) == 0


# #########  MCL stuff  ########## #
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
    assert orphans.lrg_cluster_sim_scores.iloc[0] == 0.96194166128549941
    assert orphans.lrg_cluster_sim_scores.iloc[-1] == 0.94705822908105275
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


def test_place_orphans(hf, capsys):
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
    clusters = [cluster1, cluster2, cluster3, cluster4, cluster5]

    orphans = rdmcl.Orphans(parent_sb, clusters, broker, psi_pred_ss2_dfs)

    # ##### mc_check_orphans() ##### #
    tmp_file = br.TempFile()
    orphans.mc_check_orphans(cluster3, [tmp_file.path])
    orphans.mc_check_orphans(cluster4, [tmp_file.path])
    orphans.mc_check_orphans(cluster5, [tmp_file.path])
    assert tmp_file.read() == """\
group_0_2\tgroup_0_1\t0.0335324691423
group_0_3\tgroup_0_0\t0.0147698715824
"""
    assert hf.string2hash(orphans.tmp_file.read()) == "b7c305b13ec68919ca72599cd33f1b73", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### _check_orphan() ##### #
    assert orphans._check_orphan(cluster3) == ('group_0_1', 0.033532469142303234)  # Clsuter is placed
    assert not orphans._check_orphan(cluster5)  # Insufficient support for largest cluster
    assert hf.string2hash(orphans.tmp_file.read()) == "0820336645b9080ff1f80b1d8c6fee0c", print(orphans.tmp_file.read())
    orphans.tmp_file.clear()

    # ##### place_orphans() ##### #
    orphans.place_orphans(multi_core=False)
    assert orphans.clusters[0].seq_ids == ["Hvu-PanxβA"]
    assert orphans.clusters[1].seq_ids == ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE',
                                           'Bfr-PanxαA', 'Hca-PanxαA', 'Lcr-PanxαG']
    assert orphans.clusters[2].seq_ids == ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD',
                                           'Pba-PanxαB', 'Tin-PanxαF', 'Vpa-PanxαD']
    assert hf.string2hash(orphans.tmp_file.read()) == "8c1f9952bd9f10e82688fdfde2c78b77", print(orphans.tmp_file.read())

    # Multicore doesn't seem to work in py.test, but can at least call it like it does
    orphans = rdmcl.Orphans(parent_sb, clusters, broker, psi_pred_ss2_dfs)
    orphans.place_orphans(multi_core=True)
    assert len(orphans.clusters) == 5

    # If no small clusters, nothing much happens
    orphans = rdmcl.Orphans(parent_sb, clusters[:2], broker, psi_pred_ss2_dfs)
    assert not orphans.place_orphans(multi_core=False)
    assert len(orphans.clusters) == 2

    # ToDo: Test for failure or Tukey HSD test
