#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import MyFuncs
import os
import pandas as pd
from collections import OrderedDict
from .. import rdmcl
from buddysuite import buddy_resources as br


# #########  Clusters  ########## #
def test_cluster_instantiate_group_0(hf):
    cluster = rdmcl.Cluster(hf.cteno_ids, hf.cteno_sim_scores)
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
    assert cluster.seq_ids == hf.cteno_ids
    assert cluster.seq_id_hash == hf.string2hash("".join(sorted(hf.cteno_ids)))


def test_cluster_instantiate_child(hf):
    parent = rdmcl.Cluster(hf.cteno_ids, hf.cteno_sim_scores)
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


def test_cluster_get_best_hits(hf):
    cluster = rdmcl.Cluster(hf.cteno_ids, hf.cteno_sim_scores)
    best_hit = cluster.get_best_hits("Bab-PanxαA")
    assert best_hit.iloc[0].seq2 == "Lcr-PanxαG"

# #########  PSI-PRED  ########## #
bins = ['chkparse', 'psipass2', 'psipred', 'seq2mtx']


@pytest.mark.parametrize("binary", bins)
def test_psipred_bins(binary, hf):
    assert os.path.isfile("{0}..{1}..{1}psipred{1}bin{1}{2}".format(hf.resource_path, hf.sep, binary))


def test_psi_pred(hf):
    tmpdir = br.TempDir()
    tmpdir.subdir("psi_pred")
    rdmcl.psi_pred(hf.cteno_panxs.to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}psi_pred{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.psi_pred(hf.cteno_panxs.to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "af9666d37426caa2bbf6b9075ce8df96", print(output)
