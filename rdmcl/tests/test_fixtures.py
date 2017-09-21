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
    assert str(hf.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker)) == """\
          seq1         seq2   subsmat       psi  raw_score     score
0   Lla-PanxαA   Oma-PanxαD  0.085823  0.474372   0.899475  0.202387
1   Lla-PanxαA  Mle-Panxα11  0.205013  0.369221   0.904953  0.254275
2   Lla-PanxαA   Pba-PanxαB  0.000000  0.208572   0.890717  0.062572
3  Mle-Panxα11   Oma-PanxαD  1.000000  1.000000   0.960266  1.000000
4  Mle-Panxα11   Pba-PanxαB  0.220205  0.091886   0.901932  0.181710
5   Lla-PanxαA   Tin-PanxαF  0.277117  0.262270   0.907656  0.272663
6  Mle-Panxα11   Tin-PanxαF  0.974179  0.737730   0.955062  0.903245
7   Oma-PanxαD   Pba-PanxαB  0.145473  0.191047   0.898965  0.159145
8   Pba-PanxαB   Tin-PanxαF  0.270091  0.000000   0.903550  0.189063
9   Oma-PanxαD   Tin-PanxαF  0.837363  0.729671   0.946959  0.805055"""

    assert str(hf.get_db_graph("Foo", broker)) == """\
Empty DataFrame
Columns: [seq1, seq2, subsmat, psi, raw_score, score]
Index: []"""


def test_helper_get_sim_scores(hf):
    assert str(hf.get_sim_scores(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'])) == """\
            seq1        seq2   subsmat       psi  raw_score     score
1351  Bch-PanxαC  BOL-PanxαA  0.737556  0.529124   0.727254  0.675027
5649  Bab-PanxαB  Bch-PanxαC  0.732719  0.525940   0.723317  0.670685
5666  Bab-PanxαB  BOL-PanxαA  0.956592  0.967118   0.955345  0.959750"""
