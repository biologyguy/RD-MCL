#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import pandas as pd
from collections import OrderedDict
from .. import rdmcl
from buddysuite import buddy_resources as br


# #########  Clusters  ########## #
def test_cluster_instantiate_group_0(hf):
    cluster = rdmcl.Cluster(*hf.base_cluster_args())
    assert [taxa for taxa in cluster.taxa] == ['Bab', 'Bch', 'Bfo', 'Bfr', 'BOL', 'Cfu', 'Dgl', 'Edu', 'Hca', 'Hru',
                                               'Lcr', 'Lla', 'Mle', 'Oma', 'Pba', 'Tin', 'Vpa', 'Hvu']
    assert hf.string2hash(cluster.sim_scores.to_csv()) == "984c4424c2b8529694696d715c4108a5"
    assert cluster.taxa_separator == "-"
    assert cluster.parent is None
    assert cluster.subgroup_counter == 0
    assert cluster.cluster_score is None
    assert cluster.collapsed_genes == OrderedDict([('Cfu-PanxαA', ['Cfu-PanxαF']), ('Lcr-PanxαA', ['Lcr-PanxαL']),
                                                   ('Mle-Panxα10A', ['Mle-Panxα9']),
                                                   ('Hvu-PanxβM', ['Hvu-PanxβI', 'Hvu-PanxβG', 'Hvu-PanxβH',
                                                                   'Hvu-PanxβD', 'Hvu-PanxβK', 'Hvu-PanxβE',
                                                                   'Hvu-PanxβF', 'Hvu-PanxβA', 'Hvu-PanxβC',
                                                                   'Hvu-PanxβB', 'Hvu-PanxβJ', 'Hvu-PanxβL',
                                                                   'Hvu-PanxβO'])])
    assert cluster._name == "group_0"
    assert cluster.seq_ids == ['Bab-PanxαA', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαD', 'Bab-PanxαE', 'Bch-PanxαA',
                               'Bch-PanxαB', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαA', 'Bfo-PanxαB',
                               'Bfo-PanxαC', 'Bfo-PanxαD', 'Bfo-PanxαE', 'Bfo-PanxαF', 'Bfo-PanxαG', 'Bfo-PanxαH',
                               'Bfo-PanxαI', 'Bfo-PanxαJ', 'Bfr-PanxαA', 'Bfr-PanxαB', 'Bfr-PanxαC', 'Bfr-PanxαD',
                               'BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαE', 'BOL-PanxαF',
                               'BOL-PanxαG', 'BOL-PanxαH', 'Cfu-PanxαA', 'Cfu-PanxαB', 'Cfu-PanxαC', 'Cfu-PanxαD',
                               'Cfu-PanxαE', 'Dgl-PanxαA', 'Dgl-PanxαB', 'Dgl-PanxαC', 'Dgl-PanxαD', 'Dgl-PanxαE',
                               'Dgl-PanxαF', 'Dgl-PanxαG', 'Dgl-PanxαH', 'Dgl-PanxαI', 'Edu-PanxαA', 'Edu-PanxαB',
                               'Edu-PanxαC', 'Edu-PanxαD', 'Edu-PanxαE', 'Edu-PanxαF', 'Edu-PanxαG', 'Edu-PanxαH',
                               'Edu-PanxαI', 'Hca-PanxαA', 'Hca-PanxαB', 'Hca-PanxαC', 'Hca-PanxαD', 'Hca-PanxαE',
                               'Hca-PanxαF', 'Hca-PanxαG', 'Hca-PanxαH', 'Hru-PanxαA', 'Hru-PanxαB', 'Hru-PanxαC',
                               'Hru-PanxαD', 'Hru-PanxαE', 'Lcr-PanxαA', 'Lcr-PanxαB', 'Lcr-PanxαC', 'Lcr-PanxαD',
                               'Lcr-PanxαE', 'Lcr-PanxαF', 'Lcr-PanxαG', 'Lcr-PanxαH', 'Lcr-PanxαI', 'Lcr-PanxαJ',
                               'Lcr-PanxαK', 'Lla-PanxαA', 'Lla-PanxαB', 'Lla-PanxαC', 'Mle-Panxα1', 'Mle-Panxα2',
                               'Mle-Panxα3', 'Mle-Panxα4', 'Mle-Panxα5', 'Mle-Panxα6', 'Mle-Panxα7A', 'Mle-Panxα8',
                               'Mle-Panxα10A', 'Mle-Panxα11', 'Mle-Panxα12', 'Oma-PanxαA', 'Oma-PanxαB', 'Oma-PanxαC',
                               'Oma-PanxαD', 'Pba-PanxαA', 'Pba-PanxαB', 'Pba-PanxαC', 'Pba-PanxαD', 'Pba-PanxαE',
                               'Pba-PanxαF', 'Pba-PanxαG', 'Tin-PanxαA', 'Tin-PanxαB', 'Tin-PanxαC', 'Tin-PanxαD',
                               'Tin-PanxαE', 'Tin-PanxαF', 'Vpa-PanxαA', 'Vpa-PanxαB', 'Vpa-PanxαC', 'Vpa-PanxαD',
                               'Vpa-PanxαE', 'Vpa-PanxαF', 'Vpa-PanxαG', 'Hvu-PanxβM']

    assert cluster.seq_id_hash == "1433a90160800be5348bc92477a87420"


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
    assert cluster.seq_id_hash == hf.string2hash("".join(sorted(child_ids)))


def test_cluster_get_name(hf):
    parent = rdmcl.Cluster(*hf.base_cluster_args())
    assert parent.name() == "group_0"
    child = rdmcl.Cluster(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'], parent.sim_scores, parent=parent)
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
    group_0_0 = rdmcl.Cluster(child_ids, sim_scores, parent=group_0)
    assert group_0_0._name is None

    grandchild_ids = ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB']
    group_0_0_0 = rdmcl.Cluster(grandchild_ids, sim_scores, parent=group_0_0)
    with pytest.raises(ValueError) as err:
        group_0_0_0.set_name()
    assert "Parent of current cluster has not been named." in str(err)

    group_0_0.set_name()
    assert group_0_0._name == "group_0_0"


def test_cluster_compare(hf, capsys):
    subject = rdmcl.Cluster(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE'],
                            hf.get_data("cteno_sim_scores"))
    query = rdmcl.Cluster(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE'],
                          hf.get_data("cteno_sim_scores"))
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


def test_cluster_len(hf):
    cluster = rdmcl.Cluster(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE'],
                            hf.get_data("cteno_sim_scores"))
    assert len(cluster) == 5


def test_cluster_str(hf):
    cluster = rdmcl.Cluster(['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE'],
                            hf.get_data("cteno_sim_scores"))
    assert str(cluster) == "['Hru-PanxαA', 'Lcr-PanxαH', 'Tin-PanxαC', 'Oma-PanxαC', 'Dgl-PanxαE']"

# #########  PSI-PRED  ########## #
bins = ['chkparse', 'psipass2', 'psipred', 'seq2mtx']


@pytest.mark.parametrize("binary", bins)
def test_psipred_bins(binary, hf):
    assert os.path.isfile("{0}..{1}..{1}psipred{1}bin{1}{2}".format(hf.resource_path, hf.sep, binary))


def test_psi_pred(hf):
    tmpdir = br.TempDir()
    tmpdir.subdir("psi_pred")
    rdmcl.psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}psi_pred{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.psi_pred(hf.get_data("cteno_panxs").to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "af9666d37426caa2bbf6b9075ce8df96", print(output)
