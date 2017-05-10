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
    assert hf._cteno_ids == [rec.id for rec in hf._cteno_panxs.records]
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
    assert str(hf.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker)) == """\
          seq1         seq2   subsmat       psi  raw_score     score
0  Mle-Panxα11   Oma-PanxαD  1.000000  1.000000   0.960266  1.000000
1   Lla-PanxαA   Oma-PanxαD  0.130327  0.476188   0.900593  0.234086
2  Mle-Panxα11   Tin-PanxαF  0.890459  0.660884   0.948619  0.821587
3   Lla-PanxαA   Pba-PanxαB  0.000000  0.273360   0.889770  0.082008
4   Lla-PanxαA  Mle-Panxα11  0.241332  0.395932   0.905984  0.287712
5   Pba-PanxαB   Tin-PanxαF  0.112153  0.000000   0.892309  0.078507
6   Lla-PanxαA   Tin-PanxαF  0.198788  0.208441   0.900617  0.201684
7  Mle-Panxα11   Pba-PanxαB  0.187336  0.158893   0.899186  0.178803
8   Oma-PanxαD   Pba-PanxαB  0.103258  0.285039   0.896091  0.157792
9   Oma-PanxαD   Tin-PanxαF  0.756897  0.653412   0.940558  0.725851"""

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
