#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os.path import join
import os
from .. import rename_orthogroup
from buddysuite import buddy_resources as br
from shutil import copyfile
from collections import OrderedDict
import pytest
from .. import helpers as hlp
from datetime import date


def test_prepare_clusters(hf):
    test_dir = br.TempDir()
    cluster_file = join(test_dir.path, "final_clusters.txt")
    copyfile(join(hf.resource_path, "final_clusters.txt"), cluster_file)
    with open(cluster_file, "a") as cfile:
        cfile.write("\n\n")
    expect_output = OrderedDict([('group_0_1', ['70.0358', 'BOL-PanxαC', 'Bab-PanxαD', 'Bch-PanxαE', 'Bfo-PanxαD',
                                                'Bfr-PanxαC', 'Cfu-PanxαC', 'Dgl-PanxαG', 'Edu-PanxαB', 'Hca-PanxαH']),
                                 ('group_0_3', ['45.8238', 'BOL-PanxαA', 'Bab-PanxαB', 'Bch-PanxαC', 'Bfo-PanxαB',
                                                'Dgl-PanxαE', 'Edu-PanxαA', 'Hca-PanxαB', 'Hru-PanxαA', 'Lcr-PanxαH']),
                                 ('group_0_0', ['75.9846', 'Bab-PanxαC', 'Bch-PanxαD', 'Bfo-PanxαC', 'Bfr-PanxαB',
                                                'Cfu-PanxαA', 'Cfu-PanxαD', 'Cfu-PanxαF']),
                                 ('group_0_18', ['2.0734', 'Hvu-PanxβA', 'Hvu-PanxβB', 'Hvu-PanxβC', 'Hvu-PanxβD',
                                                 'Hvu-PanxβE', 'Hvu-PanxβF', 'Hvu-PanxβG']),
                                 ('group_0_2', ['54.1235', 'BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE',
                                                'Bfr-PanxαA', 'Dgl-PanxαI', 'Edu-PanxαE']),
                                 ('group_0_5', ['20.471', 'BOL-PanxαD', 'Bab-PanxαE', 'Bfo-PanxαI', 'Dgl-PanxαD']),
                                 ('group_0_6', ['13.037', 'BOL-PanxαH', 'Dgl-PanxαH', 'Edu-PanxαC']),
                                 ('group_0_7', ['9.316', 'Bfo-PanxαG', 'Dgl-PanxαA']),
                                 ('group_0_19', ['2.9556', 'Hru-PanxαC']),
                                 ('group_0_20', ['1.642', 'Edu-PanxαI']),
                                 ('group_0_23', ['1.642', 'Edu-PanxαD']),
                                 ('group_0_26', ['2.463', 'Cfu-PanxαB']),
                                 ('group_0_30', ['1.4778', 'Bfo-PanxαA'])])

    output = rename_orthogroup.prepare_clusters(cluster_file)
    assert output == expect_output


def test_argparse_init(monkeypatch):
    rdmcl_dir = br.TempDir()
    group_name = "group_XX"
    new_name = "penguins"
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    inargs = rename_orthogroup.argparse_init()
    assert inargs.group_name == "group_XX"
    assert inargs.new_name == "penguins"
    assert inargs.force is False


def test_main(capsys, hf, monkeypatch):
    rdmcl_dir = br.TempDir()
    fake_dir = "/fake_dir"
    group_name = "group_0_1"
    fake_group_name = "group_XX"
    new_name_exists = "group_0_3"
    new_name_space = "new group"
    new_name_slash = "new/group"
    new_name = "penguins"
    new_name_too_long = "a" * 100000
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

    # Test output directory absent
    argv = ['rename_orthogroup.py', fake_dir, group_name, new_name]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert err == "Error: The provided RD-MCL output directory does not exist.\n"

    # Test final_clusters.txt file absent
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert err == "Error: The provided RD-MCL output directory does not contain the final_clusters.txt file.\n"

    # Test group_name not in clusters
    cluster_file = join(rdmcl_dir.path, "final_clusters.txt")
    copyfile(join(hf.resource_path, "final_clusters.txt"), cluster_file)
    argv = ['rename_orthogroup.py', rdmcl_dir.path, fake_group_name, new_name]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    monkeypatch.setattr(hlp, "prepare_clusters", lambda *_, **__: clusters)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "Error: The provided group name 'group_XX' does not exist." in err

    # Test new_name already exists
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name_exists]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert err == "Error: The new group name 'group_0_3' already exists, please choose a unique group name.\n"

    # Test new_name has spaces
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name_space]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "You tried to set the new name to 'new group', might I suggest 'new_group'?" in err

    # Test new_name has slashes
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name_slash]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "You tried to set the new name to 'new/group', might I suggest 'new_group'?" in err

    # Test user aborts
    argv = ['rename_orthogroup.py', rdmcl_dir.path, group_name, new_name]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    monkeypatch.setattr(br, "ask", lambda input_prompt, **__: False)
    with pytest.raises(SystemExit):
        rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "Abort!!" in out

    # Test files do not exist (except final_clusters.txt)
    monkeypatch.setattr(br, "ask", lambda input_prompt, **__: True)
    rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "/alignments/group_0_1' does not exist" in out
    assert "/hmm/group_0_1' does not exist." in out
    assert "/sim_scores/group_0_1.scores' does not exist." in out
    assert "/mcmcmc/group_0_1' does not exist." in out
    assert "/paralog_cliques' does not exist." in out
    assert "/manual_merge.log' does not exist." in out
    assert "placement.log' does not exist." in out
    assert "final_clusters.txt modified" in out

    # Create files and directories
    hmmdir = rdmcl_dir.subdir("hmm")
    hmmfile = open(join(hmmdir, group_name), 'w+').close
    aligndir = rdmcl_dir.subdir("alignments")
    alignfile = open(join(aligndir, group_name), 'w+').close
    simdir = rdmcl_dir.subdir("sim_scores")
    simfile = open(join(simdir, "{}{}".format(group_name, ".scores")), 'w+').close
    paracliques = join(rdmcl_dir.path, "paralog_cliques")
    with open(paracliques, "w") as p:
        p.write(
            """
            ###########################################################
            # If a named cluster contains reciprocal best hit cliques #
            # among a group of paralogs, they are collapsed down to a #
            # single representative. The collapses are itemized here. #
            ###########################################################
            
            # group_0_1
            {"Hvu-Panx\u03b2E": ["Hvu-Panx\u03b2M"]}
            
            ###########################################################
            #  Paralogs were expanded back into the following groups  #
            ###########################################################
            
            # group_0_0_1_1_1
            {"Hvu-PanxβE": ['Hvu-PanxβM']}
            
            """
        )
    manual_merge = join(rdmcl_dir.path, "manual_merge.log")
    with open(manual_merge, "a") as ifile:
        ifile.write("%s %s -> group_0_3\n" % (date.today(), group_name))
    placement = join(rdmcl_dir.path, "placement.log")
    with open(placement, "a") as ifile:
        ifile.write("group_0_1\n")
    rdmcl_dir.subdir("mcmcmc")
    mcdir = rdmcl_dir.subdir(join("mcmcmc", "group_0_1"))

    # Test output
    rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "/alignments/group_0_1 -->" in out
    assert os.path.isfile(join(aligndir, group_name)) is False
    assert os.path.isfile(join(aligndir, new_name)) is True
    assert "hmm/group_0_1 -->" in out
    assert os.path.isfile(join(hmmdir, group_name)) is False
    assert os.path.isfile(join(hmmdir, new_name)) is True
    assert "/sim_scores/group_0_1.scores -->" in out
    assert os.path.isfile(join(simdir, "{}{}".format(group_name, ".scores"))) is False
    assert os.path.isfile(join(simdir, "{}{}".format(new_name, ".scores"))) is True
    assert "/mcmcmc/group_0_1 -->" in out
    assert "/final_clusters.txt modified" in out
    assert "/paralog_cliques modified" in out
    assert "/manual_merge.log modified" in out
    assert "/placement.log modified" in out

    # Test file modifications
    with open(cluster_file, "r") as ifile:
        line = ifile.readline()
        assert line.startswith("penguins")
    with open(paracliques) as ifile:
        contents = ifile.read()
        assert ' # penguins' in contents
    with open(manual_merge) as ifile:
        contents = ifile.read()
        assert "penguins ->" in contents
    with open(placement) as ifile:
        contents = ifile.read()
        assert contents == "penguins\n"

    # Test name too long
    test_dir = br.TempDir()
    copyfile(join(hf.resource_path, "final_clusters.txt"), join(test_dir.path, "final_clusters.txt"))
    hmmdir = test_dir.subdir("hmm")
    hmmfile = open(join(hmmdir, group_name), 'w+').close
    aligndir = test_dir.subdir("alignments")
    alignfile = open(join(aligndir, group_name), 'w+').close
    simdir = test_dir.subdir("sim_scores")
    simfile = open(join(simdir, "{}{}".format(group_name, ".scores")), 'w+').close
    paracliques = open(join(test_dir.path, "paralog_cliques"), 'w+').close
    manual_merge = open(join(test_dir.path, "manual_merge.log"), 'w+').close
    log = open(join(test_dir.path, "placement.log"), 'w+').close
    test_dir.subdir("mcmcmc")
    mcdir = test_dir.subdir(join("mcmcmc", "group_0_1"))

    argv = ['rename_orthogroup.py', test_dir.path, group_name, new_name_too_long]
    monkeypatch.setattr(rename_orthogroup.sys, "argv", argv)
    rename_orthogroup.main()
    out, err = capsys.readouterr()
    assert "/alignments/group_0_1 -->" not in out
    assert os.path.isfile(join(aligndir, group_name)) is True
    assert os.path.isfile(join(aligndir, new_name)) is False
    assert "hmm/group_0_1 -->" not in out
    assert os.path.isfile(join(hmmdir, group_name)) is True
    assert os.path.isfile(join(hmmdir, new_name)) is False
    assert "/sim_scores/group_0_1.scores -->" not in out
    assert os.path.isfile(join(simdir, "{}{}".format(group_name, ".scores"))) is True
    assert os.path.isfile(join(simdir, "{}{}".format(new_name, ".scores"))) is False
    assert "/mcmcmc/group_0_1 -->" not in out

