#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import MyFuncs
import os
from hashlib import md5
from buddysuite import SeqBuddy as Sb
from .. import rdmcl


def string2hash(_input):
    return md5(_input.encode("utf-8")).hexdigest()

datasets = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../datasets')


# ### Begin tests
def test_psi_pred():
    seqbuddy = Sb.SeqBuddy("%s/Cteno_pannexins.fa" % datasets)
    tmpdir = MyFuncs.TempDir()
    tmpdir.subdir("psi_pred")
    rdmcl._psi_pred(seqbuddy.to_dict()["BOL-PanxαB"], [tmpdir.path])
    with open("%s/psi_pred/BOL-PanxαB.ss2" % tmpdir.path, "r") as ifile:
        assert string2hash(ifile.read()) == "4a245a9304114bc2172b865bc5a266f0"
