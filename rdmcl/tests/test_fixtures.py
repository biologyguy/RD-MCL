#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import buddysuite.SeqBuddy
from hashlib import md5


def test_helper_attributes(hf):
    assert hf.sep == os.sep
    assert hf.resource_path == os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources') + os.sep
    assert type(hf._cteno_panxs) == buddysuite.SeqBuddy.SeqBuddy
    assert hf._cteno_ids == [rec.id for rec in hf._cteno_panxs.records]
    assert md5(hf._cteno_sim_scores.to_csv().encode("utf-8")).hexdigest() == "6402f2222a0cc43d3e32c6f8cc84a47a"


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


def test_helper_get_sim_scores(hf):
    assert str(hf.get_sim_scores(['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC'])) == """\
            seq1        seq2     score
1351  Bch-PanxαC  BOL-PanxαA  0.880398
5649  Bab-PanxαB  Bch-PanxαC  0.727946
5666  Bab-PanxαB  BOL-PanxαA  0.960119"""
