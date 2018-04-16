#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from collections import OrderedDict
from os.path import join
from shutil import copyfile
from types import SimpleNamespace
from copy import copy
from datetime import date

import pandas as pd
import pytest
import scipy.stats
from buddysuite import buddy_resources as br

from .. import compare_homolog_groups
from .. import helpers as hlp
from .. import merge_orthogroups


def generate_user_answer():
    answer_list = [True, True, False, False, True, True, True, True, True, True, False]
    for answer in answer_list:
        yield answer


def mock_index_error(*args, **kwargs):
    raise IndexError(args, kwargs)


def mock_prepare_between_group_r2_df(*_, **__):
    prep_between_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαB", "Bab-PanxαB", "Bab-PanxαB"],
                            "seq2": ["BOL-PanxαA", "Dgl-PanxαE", "BOL-PanxαC", "Bab-PanxαD", "Bch-PanxαE"],
                            "r_square": [0.550504, 0.520893, 0.555286, 0.565252, 0.566482]}
    prep_between_r2_df = pd.DataFrame(data=prep_between_r2_data)
    return prep_between_r2_df


def mock_prepare_within_group_r2_df(*_, **__):
    prep_within_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαD", "Bab-PanxαD", "Bab-PanxαD"],
                           "seq2": ["Dgl-PanxαG", "Hca-PanxαH", "BOL-PanxαC", "Bch-PanxαE", "Bfo-PanxαD"],
                           "r_square": [0.997720, 0.992786, 0.995777, 0.995940, 0.995856]}
    prep_within_r2_df = pd.DataFrame(data=prep_within_r2_data)
    return prep_within_r2_df


def test_check_init(hf, monkeypatch):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    taxa = OrderedDict([('BOL', ['BOL-PanxαC', 'BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαD', 'BOL-PanxαH']),
                        ('Bab', ['Bab-PanxαD', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαA', 'Bab-PanxαE']),
                        ('Bch', ['Bch-PanxαE', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαA']),
                        ('Bfo', ['Bfo-PanxαD', 'Bfo-PanxαB', 'Bfo-PanxαC', 'Bfo-PanxαE', 'Bfo-PanxαI', 'Bfo-PanxαG',
                                 'Bfo-PanxαA']),
                        ('Bfr', ['Bfr-PanxαC', 'Bfr-PanxαB', 'Bfr-PanxαA']),
                        ('Cfu', ['Cfu-PanxαC', 'Cfu-PanxαA', 'Cfu-PanxαD', 'Cfu-PanxαF', 'Cfu-PanxαB']),
                        ('Dgl', ['Dgl-PanxαG', 'Dgl-PanxαE', 'Dgl-PanxαI', 'Dgl-PanxαD', 'Dgl-PanxαH', 'Dgl-PanxαA']),
                        ('Edu', ['Edu-PanxαB', 'Edu-PanxαA', 'Edu-PanxαE', 'Edu-PanxαC', 'Edu-PanxαI', 'Edu-PanxαD']),
                        ('Hca', ['Hca-PanxαH', 'Hca-PanxαB']),
                        ('Hru', ['Hru-PanxαA', 'Hru-PanxαC']),
                        ('Lcr', ['Lcr-PanxαH']),
                        ('Hvu', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE', 'Hvu-PanxβF',
                                 'Hvu-PanxβG'])])
    cluster_object = SimpleNamespace(taxa_separator="-", taxa=taxa)

    # First 5 lines of output from _prepare_within_group_r2_df
    prep_within_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαD", "Bab-PanxαD", "Bab-PanxαD"],
                           "seq2": ["Dgl-PanxαG", "Hca-PanxαH", "BOL-PanxαC", "Bch-PanxαE", "Bfo-PanxαD"],
                           "r_square": [0.997720, 0.992786, 0.995777, 0.995940, 0.995856]}
    prep_within_r2_df = pd.DataFrame(data=prep_within_r2_data)

    # First 5 lines of output from _prepare_within_group_fwd_df
    prep_within_fwd_data = {"hmm_id": ["Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC"],
                            "rec_id": ["Cfu-PanxαC", "Edu-PanxαB", "Bab-PanxαD", "Bch-PanxαE", "Bfo-PanxαD"],
                            "fwd_raw": [666.8759, 619.9454, 604.2560, 615.6625, 609.4811]}
    prep_within_fwd_df = pd.DataFrame(data=prep_within_fwd_data)

    # First 5 lines of output from _prepare_between_group_r2_df
    prep_between_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαB", "Bab-PanxαB", "Bab-PanxαB"],
                            "seq2": ["BOL-PanxαA", "Dgl-PanxαE", "BOL-PanxαC", "Bab-PanxαD", "Bch-PanxαE"],
                            "r_square": [0.550504, 0.520893, 0.555286, 0.565252, 0.566482]}
    prep_between_r2_df = pd.DataFrame(data=prep_between_r2_data)

    # First 5 lines of output from _prepare_between_group_fwd_df
    prep_between_fwd_data = {"hmm_id": ["Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC"],
                             "rec_id": ["Cfu-PanxαC", "Edu-PanxαB", "Bab-PanxαD", "Bch-PanxαE", "Bfo-PanxαD"],
                             "fwd_raw": [316.2845, 313.8431, 309.1125, 315.0160, 312.1405]}
    prep_between_fwd_df = pd.DataFrame(data=prep_between_fwd_data)

    # monkeypatch functions called in def __init__
    monkeypatch.setattr(hlp, "prepare_clusters", lambda *_, **__: clusters)
    monkeypatch.setattr(compare_homolog_groups, "Cluster", lambda *_: cluster_object)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_r2_df", lambda *_: prep_within_r2_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_fwd_df", lambda *_: prep_within_fwd_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_r2_df", lambda *_: prep_between_r2_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_fwd_df", lambda *_: prep_between_fwd_df)
    monkeypatch.setattr(hlp, "create_truncnorm", lambda *_: True)

    # Make tempdir. Make subdir "hmm". Copy csv files in hmm.
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "final_clusters.txt"), join(subdir, "final_clusters.txt"))
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    copyfile(join(hf.resource_path, "hmm_fwd_scores.csv"), join(subdir, "hmm_fwd_scores.csv"))

    check = merge_orthogroups.Check(rdmcl_dir=test_dir.path)
    assert check.group_name is None
    assert check.output is None
    assert check.rdmcl_dir == test_dir.path
    assert check.clusters == clusters
    assert check.master_clust.taxa == taxa
    assert check.r_squares.at[37, 'seq1'] == 'BOL-PanxαC'
    assert check.r_squares.at[316, 'seq2'] == 'Edu-PanxαC'
    assert check.r_squares['r_square'].sum() == 455.83679143893875
    assert check.fwd_scores.at[2636, 'hmm_id'] == 'Bfo-PanxαA'
    assert check.fwd_scores.at[2, 'rec_id'] == 'Edu-PanxαC'
    assert check.fwd_scores['fwd_raw'].sum() == 640373.9347999999
    assert check.within_group_r2_df.at[1, 'seq1'] == 'BOL-PanxαC'
    assert check.within_group_r2_df.at[3, 'seq2'] == 'Bch-PanxαE'
    assert check.within_group_r2_df['r_square'].sum() == 4.978079
    assert check.within_group_r2_dist is True
    assert check.within_group_fwd_df.at[1, 'hmm_id'] == 'Cfu-PanxαC'
    assert check.within_group_fwd_df.at[3, 'rec_id'] == 'Bch-PanxαE'
    assert check.within_group_fwd_df['fwd_raw'].sum() == 3116.2209
    assert type(check.within_group_fwd_dist) is scipy.stats.kde.gaussian_kde
    assert check.btw_group_r2_df.at[1, 'seq1'] == 'BOL-PanxαC'
    assert check.btw_group_r2_df.at[3, 'seq2'] == 'Bab-PanxαD'
    assert check.btw_group_r2_df['r_square'].sum() == 2.758417
    assert check.btw_group_r2_dist is True
    assert check.btw_group_fwd_df.at[1, 'hmm_id'] == 'Bch-PanxαC'
    assert check.btw_group_fwd_df.at[3, 'rec_id'] == 'Bch-PanxαE'
    assert check.btw_group_fwd_df['fwd_raw'].sum() == 1566.3966
    assert type(check.btw_group_fwd_dist) is scipy.stats.kde.gaussian_kde


# noinspection PyCallByClass
def test_prepare_within_group_r2_df(capsys, hf):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    test_dir = br.TempDir()
    test_dir.subdir("hmm")
    r_squares = pd.read_csv(join(hf.resource_path, "rsquares_matrix.csv"))

    # File "within_group_rsquares.csv" already exists
    within_group_rsquares = pd.DataFrame(columns=["seq1", "seq2", "r_square"])
    within_group_rsquares.to_csv(join(test_dir.path, "hmm", "within_group_rsquares.csv"))
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_rsquares.csv")) is True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert type(merge) is pd.DataFrame
    os.remove(join(test_dir.path, "hmm", "within_group_rsquares.csv"))

    # force=False, file "within_group_rsquares.csv" does not exist
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_rsquares.csv")) is False
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/within_group_rsquares.csv...\n"
    assert merge.at[1, 'seq1'] == 'BOL-PanxαC'
    assert merge.at[3, 'seq2'] == 'Bch-PanxαE'
    assert merge['r_square'].sum() == 124.34523271531056
    os.remove(join(test_dir.path, "hmm", "within_group_rsquares.csv"))

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=True)
    assert err == "Preparing hmm/within_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame
    assert merge.at[1, 'seq1'] == 'BOL-PanxαC'
    assert merge.at[3, 'seq2'] == 'Bch-PanxαE'
    assert merge['r_square'].sum() == 124.34523271531056
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_rsquares.csv")) is True


# noinspection PyCallByClass
def test_prepare_within_group_fwd_df(capsys, hf):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    copyfile(join(hf.resource_path, "hmm_fwd_scores.csv"), join(subdir, "hmm_fwd_scores.csv"))
    r_squares = pd.read_csv(join(test_dir.path, "hmm", "rsquares_matrix.csv"))
    fwd_scores = pd.read_csv(join(test_dir.path, "hmm", "hmm_fwd_scores.csv"))

    # force=False, file "within_group_fwd.csv" does not exist
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_fwd.csv")) is False
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/within_group_fwd.csv...\n"
    assert merge.at[1, 'hmm_id'] == 'Cfu-PanxαC'
    assert merge.at[329, 'rec_id'] == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.9229

    # File "within_group_fwd.csv" already exists, run method again
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_fwd.csv")) is True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    assert merge.at[1, 'hmm_id'] == 'Cfu-PanxαC'
    assert merge.at[329, 'rec_id'] == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.9229
    os.remove(join(test_dir.path, "hmm", "within_group_fwd.csv"))

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=True)
    assert merge.at[1, 'hmm_id'] == 'Cfu-PanxαC'
    assert merge.at[329, 'rec_id'] == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.9229


# noinspection PyCallByClass
def test_prepare_between_group_r2_df(capsys, hf):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    r_squares = pd.read_csv(join(test_dir.path, "hmm", "rsquares_matrix.csv"))

    # force=False, file "between_group_rsquares.csv" does not exist
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert merge.at[1, 'seq1'] == 'BOL-PanxαC'
    assert merge.at[978, 'seq2'] == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.82666010791604

    # Skip the main loop when file already exists
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert merge.at[1, 'seq1'] == 'BOL-PanxαC'
    assert merge.at[978, 'seq2'] == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.82666010791604

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert merge.at[1, 'seq1'] == 'BOL-PanxαC'
    assert merge.at[978, 'seq2'] == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.82666010791604


# noinspection PyCallByClass
def test_prepare_between_group_fwd_df(capsys, hf):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    copyfile(join(hf.resource_path, "hmm_fwd_scores.csv"), join(subdir, "hmm_fwd_scores.csv"))
    fwd_scores = pd.read_csv(join(test_dir.path, "hmm", "hmm_fwd_scores.csv"))

    # force=False, file "between_group_fwd.csv" does not exist
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert merge.at[1, 'hmm_id'] == 'Bch-PanxαC'
    assert merge.at[1959, 'rec_id'] == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.17520000006

    # File "between_group_fwd.csv" already exists
    assert os.path.isfile(join(test_dir.path, "hmm", "between_group_fwd.csv")) is True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert merge.at[1, 'hmm_id'] == 'Bch-PanxαC'
    assert merge.at[1959, 'rec_id'] == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.17520000006

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert merge.at[1, 'hmm_id'] == 'Bch-PanxαC'
    assert merge.at[1959, 'rec_id'] == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.17520000006


# noinspection PyCallByClass
def test_check_existing_group(hf):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    r_squares = pd.read_csv(join(test_dir.path, "hmm", "rsquares_matrix.csv"))
    group_name1 = "group_0_19"
    group_name2 = "group_X"

    # Make master_clust
    taxa = OrderedDict([('BOL', ['BOL-PanxαC', 'BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαD', 'BOL-PanxαH']),
                        ('Bab', ['Bab-PanxαD', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαA', 'Bab-PanxαE']),
                        ('Bch', ['Bch-PanxαE', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαA']),
                        ('Bfo', ['Bfo-PanxαD', 'Bfo-PanxαB', 'Bfo-PanxαC',
                                 'Bfo-PanxαE', 'Bfo-PanxαI', 'Bfo-PanxαG', 'Bfo-PanxαA']),
                        ('Bfr', ['Bfr-PanxαC', 'Bfr-PanxαB', 'Bfr-PanxαA']),
                        ('Cfu', ['Cfu-PanxαC', 'Cfu-PanxαA', 'Cfu-PanxαD', 'Cfu-PanxαF', 'Cfu-PanxαB']),
                        ('Dgl', ['Dgl-PanxαG', 'Dgl-PanxαE', 'Dgl-PanxαI', 'Dgl-PanxαD', 'Dgl-PanxαH', 'Dgl-PanxαA']),
                        ('Edu', ['Edu-PanxαB', 'Edu-PanxαA', 'Edu-PanxαE', 'Edu-PanxαC', 'Edu-PanxαI', 'Edu-PanxαD']),
                        ('Hca', ['Hca-PanxαH', 'Hca-PanxαB']),
                        ('Hru', ['Hru-PanxαA', 'Hru-PanxαC']),
                        ('Lcr', ['Lcr-PanxαH']),
                        ('Hvu', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD',
                                 'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβG'])])

    seq_ids = ['BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαC', 'BOL-PanxαD', 'BOL-PanxαH', 'Bab-PanxαA', 'Bab-PanxαB',
               'Bab-PanxαC', 'Bab-PanxαD', 'Bab-PanxαE', 'Bch-PanxαA', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαE',
               'Bfo-PanxαA', 'Bfo-PanxαB', 'Bfo-PanxαC', 'Bfo-PanxαD', 'Bfo-PanxαE', 'Bfo-PanxαG', 'Bfo-PanxαI',
               'Bfr-PanxαA', 'Bfr-PanxαB', 'Bfr-PanxαC', 'Cfu-PanxαA', 'Cfu-PanxαB', 'Cfu-PanxαC', 'Cfu-PanxαD',
               'Cfu-PanxαF', 'Dgl-PanxαA', 'Dgl-PanxαD', 'Dgl-PanxαE', 'Dgl-PanxαG', 'Dgl-PanxαH', 'Dgl-PanxαI',
               'Edu-PanxαA', 'Edu-PanxαB', 'Edu-PanxαC', 'Edu-PanxαD', 'Edu-PanxαE', 'Edu-PanxαI', 'Hca-PanxαB',
               'Hca-PanxαH', 'Hru-PanxαA', 'Hru-PanxαC', 'Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD',
               'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβG', 'Lcr-PanxαH']

    master_clust = SimpleNamespace(taxa_separator="-", max_genes_in_a_taxa=7, parent=None, seq_ids=seq_ids, taxa=taxa)

    #  Create R2 distributions
    lower = 0
    upper = 1
    within_mu = 0.857553329071
    within_sigma = 0.191460697119
    btwn_mu = 0.255164455857
    between_sigma = 0.171689140429
    within_group_r2_dist = scipy.stats.truncnorm((lower - within_mu) / within_sigma, (upper - within_mu) / within_sigma,
                                                 loc=within_mu, scale=within_sigma)
    btw_group_r2_dist = scipy.stats.truncnorm((lower - btwn_mu) / between_sigma, (upper - btwn_mu) / between_sigma,
                                              loc=btwn_mu, scale=between_sigma)
    check = SimpleNamespace(clusters=clusters, r_squares=r_squares, master_clust=master_clust,
                            within_group_r2_dist=within_group_r2_dist, btw_group_r2_dist=btw_group_r2_dist)
    merge_orthogroups.Check.check_existing_group(check, group_name=group_name1)
    assert check.output == [['group_0_7', 0.0228, 0.096, 6.319, 7.083],
                            ['group_0_18', 0.0002, 0.0907, 5.119, 5.494],
                            ['group_0_20', 0.0, 0.0, 5.056, 5.444],
                            ['group_0_23', 0.0, 0.0, 5.056, 5.444],
                            ['group_0_26', 0.0, 0.0, 5.308, 5.717],
                            ['group_0_30', 0.0, 0.0, 4.875, 5.25],
                            ['group_0_5', 0.0, 0.0018, 10.414, 11.994],
                            ['group_0_1', 0.0, 0.0025, 30.246, 34.131],
                            ['group_0_2', 0.0, 0.0041, 19.968, 22.861],
                            ['group_0_6', 0.0, 0.0191, 8.458, 9.644],
                            ['group_0_3', 0.0, 0.0193, 42.087, 38.994],
                            ['group_0_0', 0.0, 0.0312, 15.29, 17.405]]

    # Test specified group non existent
    check = SimpleNamespace(clusters=clusters, r_squares=r_squares,
                            master_clust=master_clust,
                            within_group_r2_dist=within_group_r2_dist, btw_group_r2_dist=btw_group_r2_dist)
    with pytest.raises(IndexError) as err:
        merge_orthogroups.Check.check_existing_group(check, group_name=group_name2)
    assert "Provided group name 'group_X' not found in named clusters: group_0_0" in str(err)


# noinspection PyCallByClass
def test_merge(capsys, hf, monkeypatch):
    clusters = OrderedDict([('group_0_1', ['BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD', 'Bfr-PanxαC',
                                           'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                            ('group_0_3', ['BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB', 'Dgl-PanxαE',
                                           'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                            ('group_0_0', ['Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB', 'Cfu-PanxαA',
                                           'Cfu-PanxαD', 'Cfu-PanxαF']),
                            ('group_0_18', ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE',
                                            'Hvu-PanxβF', 'Hvu-PanxβG']),
                            ('group_0_2', ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA',
                                           'Dgl-PanxαI', 'Edu-PanxαE']),
                            ('group_0_5', ['BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                            ('group_0_6', ['BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                            ('group_0_7', ['Bfo-PanxαG', 'Dgl-PanxαA']),
                            ('group_0_19', ['Hru-PanxαC']),
                            ('group_0_20', ['Edu-PanxαI']),
                            ('group_0_23', ['Edu-PanxαD']),
                            ('group_0_26', ['Cfu-PanxαB']),
                            ('group_0_30', ['Bfo-PanxαA'])])

    output = [['group_0_7', 0.022800000000000001, 0.096000000000000002, 6.319, 7.083],
              ['group_0_18', 0.00020000000000000001, 0.090700000000000003, 9.93, 7.922],
              ['group_0_20', 0.0, 0.0, 5.056, 5.444], ['group_0_23', 0.0, 0.0, 5.056, 5.444],
              ['group_0_26', 0.0, 0.0, 5.308, 5.717], ['group_0_30', 0.0, 0.0, 4.875, 5.25],
              ['group_0_5', 0.0, 0.0018, 10.414, 11.994],
              ['group_0_1', 0.0, 0.0025000000000000001, 30.246, 34.131],
              ['group_0_2', 0.0, 0.0041000000000000003, 19.968, 22.861],
              ['group_0_6', 0.0, 0.019099999999999999, 8.458, 9.644],
              ['group_0_3', 0.0, 0.019300000000000001, 42.087, 38.53],
              ['group_0_0', 0.0, 0.031199999999999999, 15.202, 17.262]]

    test_dir = br.TempDir()
    test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "final_clusters.txt"), join(test_dir.path, "final_clusters.txt"))
    open(join(test_dir.path, "hmm", 'group_0_19'), 'w+').close()
    assert os.path.isfile(join(test_dir.path, "hmm", 'group_0_19')) is True
    query_group = "group_0_19"
    merge_group_1 = "group_0_7"
    merge_group_2 = "group_0_71"
    merge_group_3 = "group_0_3"
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_r2_df", lambda *_: True)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_r2_df", lambda *_: True)
    generator = generate_user_answer()

    def monkeyask(input_prompt, **__):
        print(input_prompt)
        return next(generator)
    monkeypatch.setattr(br, "ask", monkeyask)

    # Test when user-specified merge_group is the first item in output (group_0_7) and the merge is questionable
    # User does not abort the merge
    check = SimpleNamespace(clusters=copy(clusters), output=output, group_name=query_group, rdmcl_dir=test_dir.path,
                            _prepare_within_group_r2_df=mock_prepare_within_group_r2_df,
                            _prepare_between_group_r2_df=mock_prepare_between_group_r2_df)
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_1)
    out, err = capsys.readouterr()
    assert "Less than 5% of sequences within current" in out
    log_file_path = join(test_dir.path, "manual_merge.log")
    logfile = open(log_file_path, "r")
    assert logfile.read() == "%s group_0_19 -> group_0_7\n" % date.today()
    final_clusters_file = open(join(test_dir.path, "final_clusters.txt"))
    assert final_clusters_file.read() == "group_0_1  70.0358 BOL-PanxαC Bab-PanxαD Bch-PanxαE Bfo-PanxαD Bfr-PanxαC " \
                                         "Cfu-PanxαC Dgl-PanxαG Edu-PanxαB Hca-PanxαH\n" \
                                         "group_0_3  45.8238 BOL-PanxαA Bab-PanxαB Bch-PanxαC Bfo-PanxαB Dgl-PanxαE " \
                                         "Edu-PanxαA Hca-PanxαB Hru-PanxαA Lcr-PanxαH\n" \
                                         "group_0_0  75.9846 Bab-PanxαC Bch-PanxαD Bfo-PanxαC Bfr-PanxαB Cfu-PanxαA " \
                                         "Cfu-PanxαD Cfu-PanxαF\n" \
                                         "group_0_18 2.0734  Hvu-PanxβA Hvu-PanxβB Hvu-PanxβC Hvu-PanxβD Hvu-PanxβE " \
                                         "Hvu-PanxβF Hvu-PanxβG\n" \
                                         "group_0_2  54.1235 BOL-PanxαB Bab-PanxαA Bch-PanxαA Bfo-PanxαE Bfr-PanxαA " \
                                         "Dgl-PanxαI Edu-PanxαE\n" \
                                         "group_0_5  20.471  BOL-PanxαD Bab-PanxαE Bfo-PanxαI Dgl-PanxαD\n" \
                                         "group_0_6  13.037  BOL-PanxαH Dgl-PanxαH Edu-PanxαC\n" \
                                         "group_0_7	7.083	Bfo-PanxαG	Dgl-PanxαA	Hru-PanxαC\n" \
                                         "group_0_19 2.9556  Hru-PanxαC\n" \
                                         "group_0_20 1.642   Edu-PanxαI\n" \
                                         "group_0_23 1.642   Edu-PanxαD\n" \
                                         "group_0_26 2.463   Cfu-PanxαB\n" \
                                         "group_0_30 1.4778  Bfo-PanxαA"
    final_clusters_file.close()
    assert os.path.isfile(join(test_dir.path, "hmm", 'group_0_19')) is False
    assert "Merged!" in out

    # Abort the questionable merge
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_1)
    out, err = capsys.readouterr()
    assert "Less than 5% of sequences within current" in out
    assert "Merge aborted!" in out

    # Test when user-specified merge_group is not in output (group_0_71)
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_2)
    out, err = capsys.readouterr()
    assert "Error: group_0_71 is not a group that group_0_19 can be merged with.\n" in err

    # Test when user-specified merge_group is in output, but not the first in the list (group_0_3): abort merge
    check = SimpleNamespace(clusters=copy(clusters), output=output, group_name=query_group, rdmcl_dir=test_dir.path,
                            _prepare_within_group_r2_df=mock_prepare_within_group_r2_df,
                            _prepare_between_group_r2_df=mock_prepare_between_group_r2_df)
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_3)
    out, err = capsys.readouterr()
    assert "The group that appears to be the most" in out
    assert "Merge aborted!" in out

    # Merge anyway. (User says yes to warning 1, and yes again for all other warnings
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_3)
    out, err = capsys.readouterr()
    assert "The group that appears to be the most" in out
    assert "Less than 5% of sequences within current" in out
    assert "reduce the combined orthogroup score, which means you will" in out
    assert "Last chance to abort!" in out
    logfile = open(log_file_path, "r")
    assert logfile.read() == ("%s group_0_19 -> group_0_7\n"
                              "%s group_0_19 -> group_0_3\n" % (date.today(), date.today()))
    final_clusters_file = open(join(test_dir.path, "final_clusters.txt"))
    assert final_clusters_file.read() == """\
group_0_1  70.0358 BOL-PanxαC Bab-PanxαD Bch-PanxαE Bfo-PanxαD Bfr-PanxαC Cfu-PanxαC Dgl-PanxαG Edu-PanxαB Hca-PanxαH
group_0_3	38.53	BOL-PanxαA	Bab-PanxαB	Bch-PanxαC	Bfo-PanxαB	Dgl-PanxαE	Edu-PanxαA	Hca-PanxαB	\
Hru-PanxαA	Hru-PanxαC	Lcr-PanxαH
group_0_0  75.9846 Bab-PanxαC Bch-PanxαD Bfo-PanxαC Bfr-PanxαB Cfu-PanxαA Cfu-PanxαD Cfu-PanxαF
group_0_18 2.0734  Hvu-PanxβA Hvu-PanxβB Hvu-PanxβC Hvu-PanxβD Hvu-PanxβE Hvu-PanxβF Hvu-PanxβG
group_0_2  54.1235 BOL-PanxαB Bab-PanxαA Bch-PanxαA Bfo-PanxαE Bfr-PanxαA Dgl-PanxαI Edu-PanxαE
group_0_5  20.471  BOL-PanxαD Bab-PanxαE Bfo-PanxαI Dgl-PanxαD
group_0_6  13.037  BOL-PanxαH Dgl-PanxαH Edu-PanxαC
group_0_7	7.083	Bfo-PanxαG	Dgl-PanxαA	Hru-PanxαC
group_0_19 2.9556  Hru-PanxαC
group_0_20 1.642   Edu-PanxαI
group_0_23 1.642   Edu-PanxαD
group_0_26 2.463   Cfu-PanxαB
group_0_30 1.4778  Bfo-PanxαA"""
    final_clusters_file.close()
    assert "Merged!" in out

    # Say yes at warning 1 and 2, but user changes his/her mind at warning 3
    check = SimpleNamespace(clusters=copy(clusters), output=output, group_name=query_group, rdmcl_dir=test_dir.path,
                            _prepare_within_group_r2_df=mock_prepare_within_group_r2_df,
                            _prepare_between_group_r2_df=mock_prepare_between_group_r2_df)
    merge_orthogroups.Check.merge(check, merge_group_name=merge_group_3)
    out, err = capsys.readouterr()
    assert "The group that appears to be the most" in out
    assert "Less than 5% of sequences within current" in out
    assert "reduce the combined orthogroup score, which means you will" in out
    assert "Merge aborted!" in out


# noinspection PyCallByClass
def test__str__():
    group_name = "group_0_7"
    output = [['group_0_7', 0.022800000000000001, 0.096000000000000002, 6.319, 7.083],
              ['group_0_18', 0.00020000000000000001, 0.090700000000000003, 9.93, 7.922],
              ['group_0_20', 0.0, 0.0, 5.056, 5.444], ['group_0_23', 0.0, 0.0, 5.056, 5.444],
              ['group_0_26', 0.0, 0.0, 5.308, 5.717], ['group_0_30', 0.0, 0.0, 4.875, 5.25],
              ['group_0_5', 0.0, 0.0018, 10.414, 11.994],
              ['group_0_1', 0.0, 0.0025000000000000001, 30.246, 34.131],
              ['group_0_2', 0.0, 0.0041000000000000003, 19.968, 22.861],
              ['group_0_6', 0.0, 0.019099999999999999, 8.458, 9.644],
              ['group_0_3', 0.0, 0.019300000000000001, 42.087, 38.53],
              ['group_0_0', 0.0, 0.031199999999999999, 15.202, 17.262]]

    check = SimpleNamespace(output=output, group_name=group_name)
    merge = merge_orthogroups.Check.__str__(check)
    assert merge == '\x1b[1mTesting group_0_7\x1b' \
                    '[0m\n\x1b[4mGroups      R²Within  R²Btw   OrigScore  NewScore\x1b[0m\n' \
                    'group_0_7   \x1b[91m0.0228    \x1b[91m0.096   \x1b[92m6.319      7.083\x1b[39m\n' \
                    'group_0_18  \x1b[91m0.0002    \x1b[91m0.0907  \x1b[91m9.93       7.922\x1b[39m\n' \
                    'group_0_20  \x1b[91m0.0       \x1b[92m0.0     \x1b[92m5.056      5.444\x1b[39m\n' \
                    'group_0_23  \x1b[91m0.0       \x1b[92m0.0     \x1b[92m5.056      5.444\x1b[39m\n' \
                    'group_0_26  \x1b[91m0.0       \x1b[92m0.0     \x1b[92m5.308      5.717\x1b[39m\n' \
                    'group_0_30  \x1b[91m0.0       \x1b[92m0.0     \x1b[92m4.875      5.25\x1b[39m\n' \
                    'group_0_5   \x1b[91m0.0       \x1b[92m0.0018  \x1b[92m10.414     11.994\x1b[39m\n' \
                    'group_0_1   \x1b[91m0.0       \x1b[92m0.0025  \x1b[92m30.246     34.131\x1b[39m\n' \
                    'group_0_2   \x1b[91m0.0       \x1b[92m0.0041  \x1b[92m19.968     22.861\x1b[39m\n' \
                    'group_0_6   \x1b[91m0.0       \x1b[92m0.0191  \x1b[92m8.458      9.644\x1b[39m\n' \
                    'group_0_3   \x1b[91m0.0       \x1b[92m0.0193  \x1b[91m42.087     38.53\x1b[39m\n' \
                    'group_0_0   \x1b[91m0.0       \x1b[92m0.0312  \x1b[92m15.202     17.262\x1b[39m\n'

    # test output absent
    check = SimpleNamespace(group_name=group_name, output=None)
    merge = merge_orthogroups.Check.__str__(check)
    assert str(merge) == "You must run Check.check() before printing"


def test_argparse_init(monkeypatch):
    out_dir = br.TempDir()
    group_name = "group_0_7"
    argv = ['merge_orthogroups.py', out_dir.path, group_name]
    monkeypatch.setattr(merge_orthogroups.sys, "argv", argv)
    temp_in_args = merge_orthogroups.argparse_init()
    assert temp_in_args.force is False
    assert temp_in_args.group_name == "group_0_7"
    assert temp_in_args.rdmcl_dir == out_dir.path


def test_main(capsys, hf, monkeypatch):
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "hmm_fwd_scores.csv"), join(subdir, "hmm_fwd_scores.csv"))
    group_name = "group_0_7"
    check = SimpleNamespace(group_name=None, output="", rdmcl_dir="", clusters="", master_clust="", r_squares="",
                            fwd_scores="", within_group_r2_df="", within_group_r2_dist="", within_group_fwd_df="",
                            within_group_fwd_dist="", btw_group_r2_df="", btw_group_r2_dist="", btw_group_fwd_df="",
                            btw_group_fwd_dist="", check_existing_group=lambda *_: True, merge=lambda *_: True)
    monkeypatch.setattr(merge_orthogroups, "Check", lambda *_: check)

    # Test non existing directory
    argv = ['merge_orthogroups.py', 'non/existing/dir/', group_name]
    monkeypatch.setattr(merge_orthogroups.sys, "argv", argv)
    with pytest.raises(SystemExit):
        merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert "Error: The provided RD-MCL output directory does not exist" in err

    # Test final_clusters not found
    argv = ['merge_orthogroups.py', test_dir.path, group_name]
    monkeypatch.setattr(merge_orthogroups.sys, "argv", argv)
    with pytest.raises(SystemExit):
        merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert "Error: The provided RD-MCL output directory does not contain the necessary file 'final_clusters.txt" in err

    # Test rsquares matrix not found
    copyfile(join(hf.resource_path, "final_clusters.txt"), join(test_dir.path, "final_clusters.txt"))
    with pytest.raises(SystemExit):
        merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert "Error: The provided RD-MCL output directory does not contain the necessary file 'hmm/rsquares_matrix.csv'"\
           in err

    # Test all files found, not merge
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert str(err) == ""

    # Test merge flag
    argv = ['merge_orthogroups.py', test_dir.path, group_name,  '--merge', 'group_0']
    monkeypatch.setattr(merge_orthogroups.sys, "argv", argv)
    merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert str(err) == ""

    # Test force
    argv = ['merge_orthogroups.py', test_dir.path, group_name, '--merge', 'group_0', '--force']
    monkeypatch.setattr(merge_orthogroups.sys, "argv", argv)
    merge_orthogroups.main()
    out, err = capsys.readouterr()
    assert str(err) == ""

    # Test IndexError
    check = SimpleNamespace(group_name=None, output="", rdmcl_dir="", clusters="", master_clust="", r_squares="",
                            fwd_scores="", within_group_r2_df="", within_group_r2_dist="", within_group_fwd_df="",
                            within_group_fwd_dist="", btw_group_r2_df="", btw_group_r2_dist="", btw_group_fwd_df="",
                            btw_group_fwd_dist="", check_existing_group=lambda *_: mock_index_error(),
                            merge=lambda *_: True)
    with pytest.raises(IndexError):
        merge_orthogroups.main()
