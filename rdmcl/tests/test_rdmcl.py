#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import MyFuncs
import os
from buddysuite import SeqBuddy as Sb
from .. import rdmcl


# ### Begin tests
bins = ['chkparse', 'psipass2', 'psipred', 'seq2mtx']


@pytest.mark.parametrize("binary", bins)
def test_psipred_bins(binary, hf):
    assert os.path.isfile("{0}..{1}..{1}psipred{1}bin{1}{2}".format(hf.resource_path, hf.sep, binary))


def test_psi_pred(hf):
    seqbuddy = Sb.SeqBuddy("%s/Cteno_pannexins.fa" % hf.resource_path)
    tmpdir = MyFuncs.TempDir()
    tmpdir.subdir("psi_pred")
    rdmcl.psi_pred(seqbuddy.to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("{0}{1}psi_pred{1}BOL-PanxαB.ss2".format(tmpdir.path, hf.sep), "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "da4192e125808e501495dfe4169d56c0", print(output)

    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "a") as ofile:
        ofile.write("\nfoo")

    # Confirm that sequences are not reprocessed if already present
    rdmcl.psi_pred(seqbuddy.to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "r") as ifile:
        output = ifile.read()
        assert hf.string2hash(output) == "af9666d37426caa2bbf6b9075ce8df96", print(output)
