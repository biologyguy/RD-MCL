#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from .. import compare_homolog_groups
from .. import helpers as hlp
from types import SimpleNamespace
from buddysuite import buddy_resources as br


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
