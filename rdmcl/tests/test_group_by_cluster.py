#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from .. import group_by_cluster
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br
import os
from os.path import join
import sys


def test_make_msa(hf, monkeypatch):
    seqbuddy = hf.get_data("cteno_panxs")
    seqbuddy.records = seqbuddy.records[:2]
    alb_obj = group_by_cluster.make_msa(seqbuddy, "clustalo")
    assert type(alb_obj) == Alb.AlignBuddy
    assert str(alb_obj) == """\
>Bab-PanxαA Be_abyssicola|m.8 and m.21|ML036514|937+ 2.
MLLLGSLGTIKNLSIFKDLSLDDWLDQMNRTFMFLLLCFMGTIVAVSQYTGKNISCNGFE
KFSDDFSQDYCWTQGLYTIKEAYDLPESQIPYPGIIPENVPACREHSLKNGGKIICPPPE
EIKPLTRARHLWYQWIPFYFWVIAPVFYLPYMFVKRMGLDRMKPLLKIMSDYYHCTTETP
SEEIIVKCADWVYNSIVDRLSEGSSWTSWRNRHGLGLAVLFSKLMYLGGSILVMMVTTLM
FQVGDFKTYGIEWLKQFPSDENYTTSVKHKLFPKMVACEIKRWGPSGLEEENGMCVLAPN
VIYQYIFLIMWFALAITICTNFFNIFFWVFKLTATRYTYSKLVATGHFSHKHPGWKFMYY
RIGTSGRVLLNIVAQNTNPIIFGAIMEKLTPSVIKHLRIGHVPGEYLTDPA
>Bab-PanxαB Be_abyssicola|m.19|ML47742|1063 2.
--MLDILSKFKGVTPFKGITIDDGWDQLNRSFMFVLLVVMGTTVTVRQYTGSVISCDGFK
KFGSTFAEDYCWTQGLYTVLEGYDQPSYNIPYPGLLPDELPACTPVKLKDGTRLKCPDAD
QLMSPTRISHLWYQWVPFYFWLAAAAFFMPYLLYKNFGMGDIKPLVRLLHNPVESDQ--E
LKKMTDKAATWLFYKFDLYMSEQSLVASLTRKHGLGLSMVFVKILYAAVSFCCFILTAEM
FSIGDFKTYGSKWIKKMRYEDTLATEEKDKLFPKMVACEVKRWGASGIEEEQGMCVLAPN
VINQYLFLILWFCLVFVMICNIVSIFVSLIKLLFTYGSYRRLLST-AFLRDDSAIKHMYF
NVGSSGRLILHVLANNTAPRVFEDILLTLAPKLIQRKLRGNGKAV------
"""

    seqbuddy.records = [seqbuddy.records[0]]
    alb_obj = group_by_cluster.make_msa(seqbuddy, "clustalo")
    assert type(alb_obj) == Alb.AlignBuddy
    assert str(alb_obj) == """\
>Bab-PanxαA Be_abyssicola|m.8 and m.21|ML036514|937+ 2.
MLLLGSLGTIKNLSIFKDLSLDDWLDQMNRTFMFLLLCFMGTIVAVSQYTGKNISCNGFE
KFSDDFSQDYCWTQGLYTIKEAYDLPESQIPYPGIIPENVPACREHSLKNGGKIICPPPE
EIKPLTRARHLWYQWIPFYFWVIAPVFYLPYMFVKRMGLDRMKPLLKIMSDYYHCTTETP
SEEIIVKCADWVYNSIVDRLSEGSSWTSWRNRHGLGLAVLFSKLMYLGGSILVMMVTTLM
FQVGDFKTYGIEWLKQFPSDENYTTSVKHKLFPKMVACEIKRWGPSGLEEENGMCVLAPN
VIYQYIFLIMWFALAITICTNFFNIFFWVFKLTATRYTYSKLVATGHFSHKHPGWKFMYY
RIGTSGRVLLNIVAQNTNPIIFGAIMEKLTPSVIKHLRIGHVPGEYLTDPA
"""

    # Don't modify if any sequence is reduced to nothing
    align = Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STP---YWAILP
""", in_format="fasta")

    seqbuddy = Sb.SeqBuddy(Alb.make_copy(align).records(), in_format="fasta")
    seqbuddy = Sb.clean_seq(seqbuddy)

    monkeypatch.setattr(Alb, "generate_msa", lambda *_, **__: align)
    alb_obj = group_by_cluster.make_msa(seqbuddy, "clustalo", trimal=[0.3])
    assert str(alb_obj) == str(align)

    # Don't modify if average sequence length is reduced by more than half
    align = Alb.AlignBuddy("""\
>A
MSTGTC-------
>B
M---TC-------
>C
M---TC---AILP
>D
-STPTC-YWAILP
""", in_format="fasta")

    seqbuddy = Sb.SeqBuddy(Alb.make_copy(align).records(), in_format="fasta")
    seqbuddy = Sb.clean_seq(seqbuddy)

    monkeypatch.setattr(Alb, "generate_msa", lambda *_, **__: align)
    alb_obj = group_by_cluster.make_msa(seqbuddy, "clustalo", trimal=[0.3])
    assert str(alb_obj) == str(align)

    # Remove some gaps
    alb_obj = group_by_cluster.make_msa(seqbuddy, "clustalo", trimal=[0.3, 0.55])
    assert str(alb_obj) == """\
>A
MSTGTC----
>B
M---TC----
>C
M---TCAILP
>D
-STPTCAILP
"""


def test_argparse_init(monkeypatch, hf):
    argv = ['rdmcl.py', hf.resource_path]
    monkeypatch.setattr(sys, "argv", argv)
    temp_in_args = group_by_cluster.argparse_init()
    assert temp_in_args.mode == "list"
    assert temp_in_args.aligner == "clustalo"
    assert temp_in_args.trimal == ["gappyout", 0.5, 0.75, 0.9, 0.95, "clean"]


def test_main_list(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "list"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "dffc4d7bea1917711bf01132a42f2891", print(out)


def test_sequence_file(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "--sequence_file",
            join(hf.resource_path, "Cteno_pannexins.fa")]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "dffc4d7bea1917711bf01132a42f2891", print(out)


def test_main_seqs(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "seqs"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "732ff95ca03b8d56f94318faaf746b05", print(out)


def test_main_aln(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "aln"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "b4b7dae5b4f81db7e2ea26417f85a054", print(out)


def test_main_cons(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "cons"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "564054d66d8a8927ff442dbc2918d943", print(out)


def test_main_select_groups(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "-g", "group_0_1", "group_0_3", "group_0_0"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "9a03d8619d686cc553f09138635c8685", print(out)
    assert not err

    argv = ['rdmcl.py', hf.resource_path, "-g", "group_0_1", "group_wonky"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert "'group_wonky' not present in clusters." in err


def test_main_min_max(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "-min", "5", "-max", "7"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "c15356cb75c752a1f1b4db69de58114f", print(out)


def test_main_strip_taxa(monkeypatch, hf, capsys):
    tmp_file = br.TempFile()
    seqbuddy = Sb.SeqBuddy(join(hf.resource_path, "Cteno_pannexins.fa"))
    seqbuddy = Sb.rename(seqbuddy, "^.*?\-")
    tmp_file.write(str(seqbuddy))
    argv = ['rdmcl.py', hf.resource_path, "-s", tmp_file.path, "-st"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "3020ea067affd21c77b7446f35689a6a", print(out)


def test_main_exclude_rbh_paralogs():
    pass


def test_main_include_counts(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "-ic"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "c22972e858d902a12e662e0899a3395a", print(out)

    argv = ['rdmcl.py', hf.resource_path, "cons", "-ic"]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    out, err = capsys.readouterr()
    assert hf.string2hash(out) == "48899e416bcd4272386d718de7324754", print(out)


def test_main_write(monkeypatch, hf, capsys):
    tmp_dir = br.TempDir()
    argv = ['rdmcl.py', hf.resource_path, "-w", tmp_dir.path]
    monkeypatch.setattr(sys, "argv", argv)
    group_by_cluster.main()
    root, dirs, files = next(os.walk(tmp_dir.path))
    assert sorted(files) == ['group_0_0.txt', 'group_0_1.txt', 'group_0_18.txt', 'group_0_19.txt', 'group_0_2.txt',
                             'group_0_20.txt', 'group_0_23.txt', 'group_0_26.txt', 'group_0_3.txt', 'group_0_30.txt',
                             'group_0_5.txt', 'group_0_6.txt', 'group_0_7.txt'], print(sorted(files))
    out, err = capsys.readouterr()
    assert out + err == ""


def test_main_errors(monkeypatch, hf, capsys):
    argv = ['rdmcl.py', hf.resource_path, "FOOBAR"]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        group_by_cluster.main()
    out, err = capsys.readouterr()
    assert 'Unrecognized mode, please select from ["seqs", "aln", "con", "list"].' in err

    argv = ['rdmcl.py', hf.resource_path, "aln"]
    monkeypatch.setattr(sys, "argv", argv)

    def kill1(*_):
        raise AttributeError("mock msa AttributeError fail")

    monkeypatch.setattr(group_by_cluster, "make_msa", kill1)
    with pytest.raises(SystemExit):
        group_by_cluster.main()
    out, err = capsys.readouterr()
    assert "mock msa AttributeError fail" in out

    def kill2(*_):
        raise SystemError("mock msa SystemError fail")

    monkeypatch.setattr(group_by_cluster, "make_msa", kill2)
    with pytest.raises(SystemExit):
        group_by_cluster.main()
    out, err = capsys.readouterr()
    assert "mock msa SystemError fail" in out
