#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from phylo_tools.homolog_caller import Cluster, Clusters, split_all_by_all, mcmcmc_mcl, clique_checker, orthogroup_caller, \
    merge_singles, support

@pytest.fixture(scope="session")
def cluster_sample(request):
    cluster = Cluster(["Tin-PanxαB", "Bab-PanxαE", "Bfo-PanxαI", "BOL-PanxαD", "Dgl-PanxαD", "Hru-PanxαB",
                       "Lcr-PanxαD", "Mle-Panxα2", "Pba-PanxαG"])
    return cluster


def test_init_cluster(cluster_sample):
    assert cluster_sample.cluster == ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD', 'Hru-PanxαB', 'Lcr-PanxαD',
                                      'Mle-Panxα2', 'Pba-PanxαG', 'Tin-PanxαB']
    assert cluster_sample.name == ""
    assert cluster_sample.taxa_split == "-"
    assert not cluster_sample.global_taxa_count

def test_score_cluster(cluster_sample):
    assert cluster_sample.score() == 81

