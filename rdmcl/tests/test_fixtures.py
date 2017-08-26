#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import buddysuite.SeqBuddy
from .. import helpers
from hashlib import md5
import pandas as pd


def test_helper_attributes(hf):
    assert hf.sep == os.sep
    assert hf.resource_path == os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources') + os.sep
    assert type(hf._cteno_panxs) == buddysuite.SeqBuddy.SeqBuddy
    assert hf._cteno_ids == sorted([rec.id for rec in hf._cteno_panxs.records])
    assert md5(hf._cteno_sim_scores.to_csv().encode("utf-8")).hexdigest() == "049088e80b31bac797a66534e518229f"


def test_helper_string2hash(hf):
    assert hf.string2hash("foo") == "acbd18db4cc2f85cedef654fccc4a4d8"


def test_helper_base_cluster_args(hf):
    cteno_ids, cteno_sim_scores = hf.base_cluster_args()
    assert cteno_ids == hf._cteno_ids
    assert cteno_sim_scores.to_csv() == hf._cteno_sim_scores.to_csv()


def test_helper_get_data(hf):
    assert str(hf.get_data("cteno_panxs")) == str(hf._cteno_panxs)
    assert hf.get_data("cteno_ids") == hf._cteno_ids
    assert hf.get_data("cteno_sim_scores").to_csv() == hf._cteno_sim_scores.to_csv()
    with pytest.raises(AttributeError) as err:
        hf.get_data("foo")
    assert "Unknown data type: foo" in str(err)
    assert type(hf.get_data("ss2_dfs")["Bab-PanxαA"]) == pd.DataFrame


def test_helper_get_db_graph(hf):
    broker = helpers.SQLiteBroker("%sdb.sqlite" % hf.resource_path)
    broker.start_broker()
    assert str(hf.get_db_graph("842e7bfe442b6cdc2011860514130663", broker)) == """\
         seq1        seq2   subsmat       psi  raw_score     score
0  Dgl-PanxαC  Tin-PanxαA  1.000000  1.000000   0.766242  1.000000
1  Lcr-PanxαE  Tin-PanxαA  0.486244  0.511777   0.687687  0.493904
2  Dgl-PanxαB  Tin-PanxαA  0.227795  0.446938   0.652298  0.293538
3  Dgl-PanxαC  Lcr-PanxαE  0.353235  0.476845   0.669438  0.390318
4  Dgl-PanxαB  Dgl-PanxαC  0.025460  0.331744   0.623121  0.117345
5  Dgl-PanxαB  Lcr-PanxαE  0.874808  0.520567   0.738867  0.768535
6  Bfo-PanxαH  Dgl-PanxαB  0.000000  0.000000   0.612204  0.000000
7  Bfo-PanxαH  Dgl-PanxαC  0.360949  0.474893   0.670406  0.395132
8  Bfo-PanxαH  Lcr-PanxαE  0.131930  0.274142   0.635774  0.174593
9  Bfo-PanxαH  Tin-PanxαA  0.397759  0.589353   0.677849  0.455237"""

    assert str(hf.get_db_graph("Foo", broker)) == """\
Empty DataFrame
Columns: [seq1, seq2, subsmat, psi, raw_score, score]
Index: []"""


def test_helper_get_sim_scores(hf):
    assert str(hf.get_sim_scores(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'])) == """\
            seq1        seq2   subsmat       psi  raw_score     score
1351  Bch-PanxαC  BOL-PanxαA  0.737556  0.529124   0.727254  0.675027
5649  Bab-PanxαB  Bch-PanxαC  0.732719  0.525940   0.723317  0.670685
5666  Bab-PanxαB  BOL-PanxαA  0.956592  0.967118   0.955345  0.959750""", print(str(hf.get_sim_scores(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'])))
