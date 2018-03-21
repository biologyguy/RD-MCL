#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .. import helpers as hlp
import pandas as pd
from .. import merge_orthogroups
from .. import compare_homolog_groups
import scipy.stats
from buddysuite import buddy_resources as br
from types import SimpleNamespace
from shutil import copyfile
import rdmcl
import os
from os.path import join
from buddysuite import SeqBuddy as Sb


def test_check_init(monkeypatch, hf):
    items = {1: "id1", 2: "id2"}
    d = {"r_square": [1, 2]}
    df = pd.DataFrame(data=d)
    d2 = {"hmm_id": [1, 2], "rec_id": [3, 4], "fwd_raw": [5, 6]}
    df2 = pd.DataFrame(data=d2)
    monkeypatch.setattr(hlp, "prepare_clusters", lambda *_, **__: items)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_r2_df", lambda *_: df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_within_group_fwd_df", lambda *_: df2)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_r2_df", lambda *_: df)
    monkeypatch.setattr(merge_orthogroups.Check, "_prepare_between_group_fwd_df", lambda *_: df2)
    monkeypatch.setattr(hlp, "create_truncnorm", lambda *_: True)

    #  Make tempdir. Make subdir "hmm". Copy csv files in hmm. Use that.
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    copyfile(join(hf.resource_path, "final_clusters.txt"), join(subdir, "final_clusters.txt"))
    copyfile(join(hf.resource_path, "rsquares_matrix.csv"), join(subdir, "rsquares_matrix.csv"))
    copyfile(join(hf.resource_path, "hmm_fwd_scores.csv"), join(subdir,     "hmm_fwd_scores.csv"))

    check = merge_orthogroups.Check(rdmcl_dir=test_dir.path)
    assert check.group_name is None
    assert check.output is None
    assert check.rdmcl_dir == test_dir.path
    assert check.clusters == {1: 'id1', 2: 'id2'}
    assert type(check.master_clust) == compare_homolog_groups.Cluster
    assert list(check.r_squares.columns) == ["seq1", "seq2", "r_square"]
    assert list(check.fwd_scores.columns) == ['hmm_id', 'rec_id', 'fwd_raw']
    assert check.within_group_r2_df.columns.values == 'r_square'
    assert check.within_group_r2_dist is True
    assert type(check.within_group_fwd_df) is pd.DataFrame
    assert type(check.within_group_fwd_dist) is scipy.stats.kde.gaussian_kde
    assert type(check.btw_group_r2_df) is pd.DataFrame
    assert check.btw_group_r2_dist is True
    assert type(check.btw_group_fwd_df) is pd.DataFrame
    assert type(check.btw_group_fwd_dist) is scipy.stats.kde.gaussian_kde


def test_prepare_within_group_r2_df(capsys, hf):
    clusters = {1: ["id1", "id2"], 2: ["id3", "id4"]}
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    r_squares = pd.read_csv(join(hf.resource_path, "rsquares_matrix.csv"))

    # force=False, file "within_group_rsquares.csv" does not exist, seqs > 1
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/within_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame

    # Create file, run method again
    within_group_rsquares = pd.DataFrame(columns=["seq1", "seq2", "r_square"])
    within_group_rsquares.to_csv(join(test_dir.path, "hmm", "within_group_rsquares.csv"))
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert type(merge) is pd.DataFrame

    # force=True, seqs < 2
    clusters2 = {1: ["id1"]}
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters2, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_within_group_r2_df(check, force=True)
    assert type(merge) is pd.DataFrame


def test_prepare_within_group_fwd_df(capsys, hf):
    clusters = {1: ["id1", "id2"], 2: ["id3", "id4"]}
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    r_squares = pd.read_csv(join(hf.resource_path, "rsquares_matrix.csv"))
    data = {"hmm_id": [1, 2], "rec_id": [3, 4], "fwd_raw": [5, 6]}
    fwd_scores = pd.DataFrame(data=data)

    # force=False, file "within_group_fwd.csv" does not exist, seqs > 1
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/within_group_fwd.csv...\n"
    assert type(merge) is pd.DataFrame

    # Create file "within_group_fwd.csv", run method again
    fwd_scores.to_csv(join(test_dir.path, "hmm", "within_group_fwd.csv"))
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert type(merge) is pd.DataFrame

    # force=True, seqs < 2
    clusters2 = {1: ["id1"]}
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters2, r_squares=r_squares, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_within_group_fwd_df(check, force=True)
    out, err = capsys.readouterr()
    assert type(merge) is pd.DataFrame


def test_prepare_between_group_r2_df(hf, capsys):
    clusters = {1: ["id1", "id2"], 2: ["id3", "id4"], 3: ["id5", "id6", "id7"]}
    clusters2 = {1: ["id1"]}
    clusters3 = {1: ["id1", "id2"], 2: ["id3"]}

    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    r_squares = pd.read_csv(join(hf.resource_path, "rsquares_matrix.csv"))

    # force=False, file "between_group_rsquares.csv" does not exist, seqs1 > 1, seqs2 >1
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame

    # force=True, seqs1 < 2 seqs2 >1
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters2, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame

    # Need to enter: if len(seqs2) < 2: in second for loop
    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters3, r_squares=r_squares)
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_rsquares.csv...\n"
    assert type(merge) is pd.DataFrame

    # Skip the main loop when file already exists
    merge = merge_orthogroups.Check._prepare_between_group_r2_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert type(merge) is pd.DataFrame


def test_prepare_between_group_fwd_df(hf, capsys):
    test_dir = br.TempDir()
    subdir = test_dir.subdir("hmm")
    clusters = {1: ["id1", "id2"], 2: ["id3", "id4"], 3: ["id5", "id6", "id7"]}
    clusters2 = {1: ["id1"]}
    clusters3 = {1: ["id1", "id2"], 2: ["id3"]}
    data = {"hmm_id": [1, 2], "rec_id": [3, 4], "fwd_raw": [5, 6]}
    fwd_scores = pd.DataFrame(data=data)

    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert type(merge) is pd.DataFrame

    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters2, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert type(merge) is pd.DataFrame

    check = SimpleNamespace(rdmcl_dir=test_dir.path, clusters=clusters3, fwd_scores=fwd_scores)
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=True)
    out, err = capsys.readouterr()
    assert err == "Preparing hmm/between_group_fwd.csv...\n"
    assert type(merge) is pd.DataFrame

    # Skip the main loop when file already exists
    merge = merge_orthogroups.Check._prepare_between_group_fwd_df(check, force=False)
    out, err = capsys.readouterr()
    assert err == ""
    assert type(merge) is pd.DataFrame



###########  FIX ISSUE WITH buddy_resources  #############

#def test_check_existing_group(hf, capsys):
#    clusters = {1: ["id1", "id2"], 2: ["id3", "id4"], 3: ["id5", "id6", "id7"]}
#    r_squares = pd.read_csv(join(hf.resource_path, "rsquares_matrix.csv"))
#    group_name = "group_0"

#    check = SimpleNamespace(clusters=clusters, group_name=group_name, r_squares=r_squares)
#    merge = merge_orthogroups.Check.check_existing_group(check, group_name=group_name)




def test_mc_fwd_back_old_hmms(monkeypatch, capsys, hf):
    seqs_file = join(hf.resource_path, "input_seqs.fa")
    sequences = Sb.SeqBuddy(seqs_file)


    query_file = br.TempFile()
    query_file.write(rec.format("fasta"))
    sequences = Sb.SeqBuddy(seqs_file)
    seq_chuncks = hlp.chunk_list(sequences.records, br.usable_cpu_count())


    seq1 = Sb.
    seq_chunk = [seq1, seq2, seq3, seq4]




    hmm_scores_file = "a"
    hmm_dir_path = "b"
    query_file = "c"
    y = [hmm_scores_file, hmm_dir_path, query_file]

    check = merge_orthogroups.Check(rdmcl_dir=test_dir.path)
    merge = merge_orthogroups.Check._mc_fwd_back_old_hmms(seq_chunk, args)
























