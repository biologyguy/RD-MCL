#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import buddysuite.SeqBuddy
from rdmcl import helpers
from hashlib import md5
import pandas as pd

pd.set_option('expand_frame_repr', False)


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
    graph = hf.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker)
    assert str(graph) == """\
          seq1         seq2         subsmat             psi       raw_score           score
0   Lla-PanxαA   Oma-PanxαD  0.085822632066  0.474372002639  0.899474740758  0.202387443238
1   Lla-PanxαA  Mle-Panxα11  0.205012941498  0.369221443725  0.904952928624  0.254275492166
2   Lla-PanxαA   Pba-PanxαB  0.000000000000  0.208571785135  0.890717163169  0.062571535540
3  Mle-Panxα11   Oma-PanxαD  1.000000000000  1.000000000000  0.960265758290  1.000000000000
4  Mle-Panxα11   Pba-PanxαB  0.220205363308  0.091886398888  0.901931544743  0.181709673982
5   Lla-PanxαA   Tin-PanxαF  0.277116813973  0.262270239432  0.907656135299  0.272662841611
6  Mle-Panxα11   Tin-PanxαF  0.974179389024  0.737729760568  0.955061737042  0.903244500487
7   Oma-PanxαD   Pba-PanxαB  0.145473499447  0.191046691983  0.898965077972  0.159145457208
8   Pba-PanxαB   Tin-PanxαF  0.270090590010  0.000000000000  0.903549609015  0.189063413007
9   Oma-PanxαD   Tin-PanxαF  0.837362886493  0.729671426789  0.946958793306  0.805055448582""", print(graph)

    assert str(hf.get_db_graph("Foo", broker)) == """\
Empty DataFrame
Columns: [seq1, seq2, subsmat, psi, raw_score, score]
Index: []"""


def test_helper_get_sim_scores(hf):
    sim_scores = hf.get_sim_scores(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'])
    assert str(sim_scores) == """\
            seq1        seq2         subsmat             psi       raw_score           score
1351  Bch-PanxαC  BOL-PanxαA  0.737556323241  0.529123922436  0.727253818704  0.675026602999
5649  Bab-PanxαB  Bch-PanxαC  0.732718501569  0.525940464746  0.723317290588  0.670685090522
5666  Bab-PanxαB  BOL-PanxαA  0.956592383247  0.967117942419  0.955344968094  0.959750050998""", print(sim_scores)
