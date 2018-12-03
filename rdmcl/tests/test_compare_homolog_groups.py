#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from rdmcl import compare_homolog_groups
from rdmcl import helpers as hlp
from types import SimpleNamespace
from buddysuite import buddy_resources as br
import sys


def test_comparison_init(monkeypatch, capsys):
    monkeypatch.setattr(compare_homolog_groups.hlp, "prepare_clusters",
                        lambda ifile: [[print("prepare_clusters => %s" % ifile), "1", "2", "3"], [1, 2, 3]])
    monkeypatch.setattr(compare_homolog_groups.Comparison, "_prepare_difference",
                        lambda *_: print("_prepare_difference()"))
    comp = compare_homolog_groups.Comparison("file_1", "file_2")
    out, err = capsys.readouterr()
    assert "prepare_clusters => file_1" in out
    assert "prepare_clusters => file_2" in out
    assert "_prepare_difference()" in out

    assert comp.true_clusters
    assert comp.query_clusters
    assert comp.total_size == 7
    assert comp.precision == 0
    assert comp.recall == 0
    assert comp.accuracy == 0
    assert comp.tn_rate == 0
    assert comp.query_score == 0
    assert comp.true_score == 0
    assert comp.pretty_out == ""


def test_comarison_prepare_difference():
    true = br.TempFile()
    true.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Bfo-PanxαB
BOL-PanxαH Dgl-PanxαH Edu-PanxαC
Bfo-PanxαG Dgl-PanxαA
Hru-PanxαC
""")
    true_clusters = hlp.prepare_clusters(true.path)

    query = br.TempFile()
    query.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Bfo-PanxαB Edu-PanxαC
BOL-PanxαH Dgl-PanxαH
Bfo-PanxαG Dgl-PanxαA
Hru-PanxαC""")
    query_clusters = hlp.prepare_clusters(query.path)

    comp = SimpleNamespace(true_clusters=true_clusters, query_clusters=query_clusters,
                           pretty_out="", _metrics=lambda *_: True,
                           _prepare_difference=compare_homolog_groups.Comparison._prepare_difference)
    comp._prepare_difference(comp)
    assert comp.pretty_out == """\
{0}BOL-PanxαA{1}\t{0}Bab-PanxαB{1}\t{0}Bch-PanxαC{1}\t{0}Bfo-PanxαB{1}\t{2}Edu-PanxαC{1}
{0}BOL-PanxαH{1}\t{0}Dgl-PanxαH{1}
{0}Bfo-PanxαG{1}\t{0}Dgl-PanxαA{1}
{0}Hru-PanxαC{1}
""".format(hlp.GREEN, hlp.DEF_FONT, hlp.RED), print(comp.pretty_out)

    true.clear()
    true.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Mle-Panxα6
Bfo-PanxαG Dgl-PanxαA Mle-Panxα12
Hru-PanxαC
Mle-Panxα1
""")

    query = br.TempFile()
    query.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Mle-Panxα6 Dgl-PanxαA
Bfo-PanxαG Mle-Panxα12
Hru-PanxαC
Mle-Panxα1""")

    comp.true_clusters = hlp.prepare_clusters(true.path)
    comp.query_clusters = hlp.prepare_clusters(query.path)
    comp.pretty_out = ""
    comp._prepare_difference(comp)

    assert comp.pretty_out == """\
{0}BOL-PanxαA{1}\t{0}Bab-PanxαB{1}\t{0}Bch-PanxαC{1}\t{0}{2}Mle-Panxα6\033[24m{1}\t{3}Dgl-PanxαA{1}
{0}Bfo-PanxαG{1}\t{0}{2}Mle-Panxα12\033[24m{1}
{0}Hru-PanxαC{1}
{0}{2}Mle-Panxα1\033[24m{1}
""".format(hlp.GREEN, hlp.DEF_FONT, hlp.UNDERLINE, hlp.RED), print(comp.pretty_out)

    query.clear()
    query.write("""\
Bab-PanxαB Bch-PanxαC Bfo-PanxαG
Mle-Panxα6 BOL-PanxαA
Hru-PanxαC Mle-Panxα1
Dgl-PanxαA
Mle-Panxα12""")

    comp.query_clusters = hlp.prepare_clusters(query.path)
    comp.pretty_out = ""
    comp._prepare_difference(comp)

    assert comp.pretty_out == """\
{3}Bab-PanxαB{1}\t{3}Bch-PanxαC{1}\t{3}Bfo-PanxαG{1}
{0}{2}Mle-Panxα6\033[24m{1}\t{0}BOL-PanxαA{1}
{0}Hru-PanxαC{1}\t{3}{2}Mle-Panxα1\033[24m{1}
{0}Dgl-PanxαA{1}
{3}{2}Mle-Panxα12\033[24m{1}
""".format(hlp.GREEN, hlp.DEF_FONT, hlp.UNDERLINE, hlp.RED), print(comp.pretty_out)

    with pytest.raises(ValueError) as err:
        query.clear()
        query.write("""\
    Bab-PanxαB Bch-PanxαC Bfo-PanxαG
    Mle-Panxα6 BOL-PanxαA
    Hru-PanxαC Mle-Panxα1""")

        comp.query_clusters = hlp.prepare_clusters(query.path)
        comp._prepare_difference(comp)
    assert "Not all sequences are present in both sets of clusters" in str(err)


def test_comparison_metrics():
    true = br.TempFile()
    true.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Bfo-PanxαB
BOL-PanxαH Dgl-PanxαH Edu-PanxαC
Bfo-PanxαG Dgl-PanxαA
Hru-PanxαC
""")
    true_clusters = hlp.prepare_clusters(true.path)

    query = br.TempFile()
    query.write("""\
BOL-PanxαA Bab-PanxαB Bch-PanxαC Bfo-PanxαB Edu-PanxαC
BOL-PanxαH Dgl-PanxαH
Bfo-PanxαG Dgl-PanxαA
Hru-PanxαC""")
    query_clusters = hlp.prepare_clusters(query.path)

    comp = SimpleNamespace(true_clusters=true_clusters, query_clusters=query_clusters, total_size=10,
                           precision=0, recall=0, accuracy=0, tn_rate=0, query_score=0, true_score=0,
                           pretty_out="", _metrics=compare_homolog_groups.Comparison._metrics)
    true_to_query = [[0, 0], [1, 1], [2, 2], [3, 3]]
    comp._metrics(comp, true_to_query)
    assert comp.precision == 0.9
    assert comp.recall == 0.9
    assert comp.accuracy == 0.95
    assert round(comp.tn_rate, 3) == 0.967

    assert comp.query_score > comp.true_score
    assert round(comp.query_score, 2) == 17.46
    assert round(comp.true_score, 2) == 16.65


def test_comparison_score():
    comp = SimpleNamespace(precision=0.1, recall=0.2, accuracy=0.3, tn_rate=0.4, query_score=100, true_score=200,
                           score=compare_homolog_groups.Comparison.score)
    score = comp.score(comp)
    assert score == """\
Precision:    10.0%
Recall:       20.0%
Accuracy:     30.0%
tn rate:      40.0%
Query score:  100
True score:   200
""", print(score)


def test_comparison_str():
    comp = SimpleNamespace(precision=0.1, recall=0.2, accuracy=0.3, tn_rate=0.4, query_score=100, true_score=200,
                           score=lambda *_: "Printing out a Comparison object.",
                           __str__=compare_homolog_groups.Comparison.__str__)

    assert comp.__str__(comp) == "Printing out a Comparison object."


def test_main(monkeypatch, capsys):
    class MockComp(object):
        def __init__(self, true_clusts, false_clusts):
            assert true_clusts
            assert false_clusts
            self.pretty_out = "Pretty Out!"

        @staticmethod
        def score():
            return "Mock score"

    argv = ['compare_homolog_groups.py', "true_clust", "false_clust"]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(compare_homolog_groups, "Comparison", MockComp)
    compare_homolog_groups.main()
    out, err = capsys.readouterr()

    assert "Pretty Out!" in out

    argv = ['compare_homolog_groups.py', "true_clust", "false_clust", "-s"]
    monkeypatch.setattr(sys, "argv", argv)
    monkeypatch.setattr(compare_homolog_groups, "Comparison", MockComp)
    compare_homolog_groups.main()
    out, err = capsys.readouterr()

    assert "Mock score" in out


def test_main_errors(monkeypatch, capsys):
    def mock_raise_valueerror1(*_):
        raise ValueError("Not all sequences are present")

    def mock_raise_valueerror2(*_):
        raise ValueError("Foo Bar Baz")

    argv = ['compare_homolog_groups.py', "true_clust", "false_clust"]
    monkeypatch.setattr(sys, "argv", argv)

    monkeypatch.setattr(compare_homolog_groups, "Comparison", mock_raise_valueerror1)
    compare_homolog_groups.main()
    out, err = capsys.readouterr()

    assert "There are differences in the sequences present in your query" in err

    monkeypatch.setattr(compare_homolog_groups, "Comparison", mock_raise_valueerror2)
    with pytest.raises(ValueError) as err:
        compare_homolog_groups.main()

    assert "Foo Bar Baz" in str(err)
