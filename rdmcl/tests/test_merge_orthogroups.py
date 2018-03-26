#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .. import helpers as hlp
import pandas as pd
from .. import merge_orthogroups
from .. import compare_homolog_groups
import numpy as np
from scipy import stats
from buddysuite import buddy_resources as br
from types import SimpleNamespace
from shutil import copyfile
import rdmcl
import os
from os.path import join
from buddysuite import SeqBuddy as Sb
from collections import OrderedDict


def test_check_init(monkeypatch, hf):
    #  Clusters corresponding to final_clusters.txt in unit_test_resources
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

    #  Make SimpleNameSpace withs some attributes of a cluster object
    taxa = OrderedDict([('BOL', ['BOL-PanxαC', 'BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαD', 'BOL-PanxαH']),
                        ('Bab', ['Bab-PanxαD', 'Bab-PanxαB', 'Bab-PanxαC','Bab-PanxαA', 'Bab-PanxαE']),
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

    #  First 5 lines of output from _prepare_within_group_r2_df
    prep_within_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαD", "Bab-PanxαD", "Bab-PanxαD"],
                           "seq2": ["Dgl-PanxαG", "Hca-PanxαH", "BOL-PanxαC", "Bch-PanxαE", "Bfo-PanxαD"],
                           "r_square": [0.997720, 0.992786, 0.995777, 0.995940, 0.995856]}
    prep_within_r2_df = pd.DataFrame(data=prep_within_r2_data)

    #  First 5 lines of output from _prepare_within_group_fwd_df
    prep_within_fwd_data = {"hmm_id": ["Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC", "Cfu-PanxαC"],
                            "rec_id": ["Cfu-PanxαC", "Edu-PanxαB", "Bab-PanxαD", "Bch-PanxαE", "Bfo-PanxαD"],
                            "fwd_raw": [666.8759, 619.9454, 604.2560, 615.6625, 609.4811]}
    prep_within_fwd_df = pd.DataFrame(data=prep_within_fwd_data)

    #  First 5 lines of output from _prepare_between_group_r2_df
    prep_between_r2_data = {"seq1": ["BOL-PanxαC", "BOL-PanxαC", "Bab-PanxαB", "Bab-PanxαB", "Bab-PanxαB"],
                            "seq2": ["BOL-PanxαA", "Dgl-PanxαE", "BOL-PanxαC", "Bab-PanxαD", "Bch-PanxαE"],
                            "r_square": [0.550504, 0.520893, 0.555286, 0.565252, 0.566482]}
    prep_between_r2_df = pd.DataFrame(data=prep_between_r2_data)

    #  First 5 lines of output from _prepare_between_group_fwd_df
    prep_between_fwd_data = {"hmm_id": ["Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC", "Bch-PanxαC"],
                             "rec_id": ["Cfu-PanxαC", "Edu-PanxαB", "Bab-PanxαD", "Bch-PanxαE", "Bfo-PanxαD"],
                             "fwd_raw": [316.2845, 313.8431, 309.1125, 315.0160, 312.1405]}
    prep_between_fwd_df = pd.DataFrame(data=prep_between_fwd_data)

    #  monkeypatch functions called in def __init__
    monkeypatch.setattr(hlp, "prepare_clusters", lambda *_, **__: clusters)
    monkeypatch.setattr(compare_homolog_groups, "Cluster", lambda *_: cluster_object)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_r2_df", lambda *_: prep_within_r2_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_fwd_df", lambda *_: prep_within_fwd_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_r2_df", lambda *_: prep_between_r2_df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_fwd_df", lambda *_: prep_between_fwd_df)
    monkeypatch.setattr(hlp, "create_truncnorm", lambda *_: True)

    #  Make tempdir. Make subdir "hmm". Copy csv files in hmm.
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
    assert check.r_squares.get_value(37, 'seq1') == 'BOL-PanxαC'
    assert check.r_squares.get_value(316, 'seq2') == 'Edu-PanxαC'
    assert check.r_squares['r_square'].sum() == 455.83679143893914
    assert check.fwd_scores.get_value(2636, 'hmm_id') == 'Bfo-PanxαA'
    assert check.fwd_scores.get_value(2, 'rec_id') == 'Edu-PanxαC'
    assert check.fwd_scores['fwd_raw'].sum() == 640373.9347999988
    assert check.within_group_r2_df.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert check.within_group_r2_df.get_value(3, 'seq2') == 'Bch-PanxαE'
    assert check.within_group_r2_df['r_square'].sum() == 4.978079
    assert check.within_group_r2_dist is True
    assert check.within_group_fwd_df.get_value(1, 'hmm_id') == 'Cfu-PanxαC'
    assert check.within_group_fwd_df.get_value(3, 'rec_id') == 'Bch-PanxαE'
    assert check.within_group_fwd_df['fwd_raw'].sum() == 3116.2209
    assert type(check.within_group_fwd_dist) is scipy.stats.kde.gaussian_kde
    assert check.btw_group_r2_df.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert check.btw_group_r2_df.get_value(3, 'seq2') == 'Bab-PanxαD'
    assert check.btw_group_r2_df['r_square'].sum() == 2.758417
    assert check.btw_group_r2_dist is True
    assert check.btw_group_fwd_df.get_value(1, 'hmm_id') == 'Bch-PanxαC'
    assert check.btw_group_fwd_df.get_value(3, 'rec_id') == 'Bch-PanxαE'
    assert check.btw_group_fwd_df['fwd_raw'].sum() == 1566.3966
    assert type(check.btw_group_fwd_dist) is scipy.stats.kde.gaussian_kde


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
    subdir = test_dir.subdir("hmm")
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
    assert merge.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert merge.get_value(3, 'seq2') == 'Bch-PanxαE'
    assert merge['r_square'].sum() == 124.34523271531056
    os.remove(join(test_dir.path, "hmm", "within_group_rsquares.csv"))

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=True)
    assert err == "Preparing hmm/within_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame
    assert merge.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert merge.get_value(3, 'seq2') == 'Bch-PanxαE'
    assert merge['r_square'].sum() == 124.34523271531056
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_rsquares.csv")) is True


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

    #  force=False, file "within_group_fwd.csv" does not exist
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_fwd.csv")) is False
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/within_group_fwd.csv...\n"
    assert merge.get_value(1, 'hmm_id') == 'Cfu-PanxαC'
    assert merge.get_value(329, 'rec_id') == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.92289999998

    #  File "within_group_fwd.csv" already exists, run method again
    assert os.path.isfile(join(test_dir.path, "hmm", "within_group_fwd.csv")) is True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    assert merge.get_value(1, 'hmm_id') == 'Cfu-PanxαC'
    assert merge.get_value(329, 'rec_id') == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.92289999998
    os.remove(join(test_dir.path, "hmm", "within_group_fwd.csv"))

    #  force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=True)
    assert merge.get_value(1, 'hmm_id') == 'Cfu-PanxαC'
    assert merge.get_value(329, 'rec_id') == 'BOL-PanxαH'
    assert merge['fwd_raw'].sum() == 177372.92289999998


def test_prepare_between_group_r2_df(hf, capsys):
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

    #  force=False, file "between_group_rsquares.csv" does not exist
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert merge.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert merge.get_value(978, 'seq2') == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.8266601079159

    #  Skip the main loop when file already exists
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert merge.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert merge.get_value(978, 'seq2') == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.8266601079159

    # force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert merge.get_value(1, 'seq1') == 'BOL-PanxαC'
    assert merge.get_value(978, 'seq2') == 'Dgl-PanxαA'
    assert merge['r_square'].sum() == 250.8266601079159


def test_prepare_between_group_fwd_df(hf, capsys):
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

    #  force=False, file "between_group_fwd.csv" does not exist
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert merge.get_value(1, 'hmm_id') == 'Bch-PanxαC'
    assert merge.get_value(1959, 'rec_id') == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.1751999998

    #  File "between_group_fwd.csv" already exists
    assert os.path.isfile(join(test_dir.path, "hmm", "between_group_fwd.csv")) is True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert merge.get_value(1, 'hmm_id') == 'Bch-PanxαC'
    assert merge.get_value(1959, 'rec_id') == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.1751999998

    #  force=True
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert merge.get_value(1, 'hmm_id') == 'Bch-PanxαC'
    assert merge.get_value(1959, 'rec_id') == 'Bfo-PanxαG'
    assert merge['fwd_raw'].sum() == 384849.1751999998


def test_check_existing_group(hf, capsys):
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

    # Make master_clust
    taxa = OrderedDict([('BOL', ['BOL-PanxαC', 'BOL-PanxαA', 'BOL-PanxαB', 'BOL-PanxαD', 'BOL-PanxαH']),
                 ('Bab', ['Bab-PanxαD', 'Bab-PanxαB', 'Bab-PanxαC', 'Bab-PanxαA', 'Bab-PanxαE']),
                 ('Bch', ['Bch-PanxαE', 'Bch-PanxαC', 'Bch-PanxαD', 'Bch-PanxαA']), ('Bfo', ['Bfo-PanxαD', 'Bfo-PanxαB',
                                                                                             'Bfo-PanxαC', 'Bfo-PanxαE',
                                                                                             'Bfo-PanxαI', 'Bfo-PanxαG',
                                                                                             'Bfo-PanxαA']),
                 ('Bfr', ['Bfr-PanxαC', 'Bfr-PanxαB', 'Bfr-PanxαA']),
                 ('Cfu', ['Cfu-PanxαC', 'Cfu-PanxαA', 'Cfu-PanxαD', 'Cfu-PanxαF', 'Cfu-PanxαB']),
                 ('Dgl', ['Dgl-PanxαG', 'Dgl-PanxαE', 'Dgl-PanxαI', 'Dgl-PanxαD', 'Dgl-PanxαH', 'Dgl-PanxαA']),
                 ('Edu', ['Edu-PanxαB', 'Edu-PanxαA', 'Edu-PanxαE', 'Edu-PanxαC', 'Edu-PanxαI', 'Edu-PanxαD']),
                 ('Hca', ['Hca-PanxαH', 'Hca-PanxαB']), ('Hru', ['Hru-PanxαA', 'Hru-PanxαC']), ('Lcr', ['Lcr-PanxαH']),
                 ('Hvu',
                  ['Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD', 'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβG'])])

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
    between_mu = 0.255164455857
    between_sigma = 0.171689140429
    within_group_r2_dist = stats.truncnorm((lower - within_mu) / within_sigma, (upper - within_mu) / within_sigma,
                                           loc=within_mu, scale=within_sigma)
    btw_group_r2_dist = stats.truncnorm((lower - between_mu) / between_sigma, (upper - between_mu) / between_sigma,
                                        loc=between_mu, scale=between_sigma)
    check = SimpleNamespace(clusters=clusters, group_name=group_name1, r_squares=r_squares, master_clust=master_clust,
                            within_group_r2_dist=within_group_r2_dist, btw_group_r2_dist=btw_group_r2_dist)
    merge = merge_orthogroups.Check.check_existing_group(check, group_name=group_name1)
    out, err = capsys.readouterr()
    assert check.output == [['group_0_7', 0.022800000000000001, 0.096000000000000002, 6.319, 7.083],
                            ['group_0_18', 0.00020000000000000001, 0.090700000000000003, 9.93, 7.922],
                            ['group_0_20', 0.0, 0.0, 5.056, 5.444], ['group_0_23', 0.0, 0.0, 5.056, 5.444],
                            ['group_0_26', 0.0, 0.0, 5.308, 5.717], ['group_0_30', 0.0, 0.0, 4.875, 5.25],
                            ['group_0_5', 0.0, 0.0018, 10.414, 11.994],
                            ['group_0_1', 0.0, 0.0025000000000000001, 30.246, 34.131],
                            ['group_0_2', 0.0, 0.0041000000000000003, 19.968, 22.861],
                            ['group_0_6', 0.0, 0.019099999999999999, 8.458, 9.644],
                            ['group_0_3', 0.0, 0.019300000000000001, 42.087, 38.53],
                            ['group_0_0', 0.0, 0.031199999999999999, 15.202, 17.262]]

    










