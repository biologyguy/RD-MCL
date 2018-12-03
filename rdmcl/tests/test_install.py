#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from rdmcl import install
import shutil
from buddysuite import buddy_resources as br
import urllib.request
import os
from os.path import join


class MockPopen(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    @staticmethod
    def wait():
        return


def test_do_not_install(capsys, monkeypatch):
    # In this case, it doesn't matter if SCRIPT_PATH is monkeypatched or not
    cwd = os.getcwd()
    temp_dir = br.TempDir()
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)

    # Monkeypatch which to show nothing is installed. User does not want to install
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)
    monkeypatch.setattr(shutil, "which", lambda *_, **__: False)
    monkeypatch.setattr(br, "ask", lambda *_, **__: False)

    install.setup()
    out, err = capsys.readouterr()
    assert "\033[1mChecking for PSIPRED:\033[m \033[91mMissing\033[39m\n\n"
    assert "\033[91mRD-MCL depends on PSIPRED, and it is not installed correctly.\033[39m\n" in out
    os.chdir(cwd)


def test_fail_psipred_local_install(capsys, monkeypatch):
    # Which returns false
    # User wants to install psipred
    # Conda is not found
    # Downloads psipred from URL
    # Installation fails for some reason, program exits
    cwd = os.getcwd()
    temp_dir = br.TempDir()

    monkeypatch.setattr(br, "TempDir", lambda *_: temp_dir)
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)

    # Nothing is installed. User agrees to install.
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)
    monkeypatch.setattr(shutil, "which", lambda *_, **__: False)
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    # Pretend to download psipred-4.01-1.tar.bz2
    monkeypatch.setattr(urllib.request, "urlretrieve", lambda *_, **__: True)

    # Create required directories for local install
    temp_dir.subdir("bin")
    temp_dir.subdir("psipred-4.01-1")
    temp_dir.subdir("share")
    temp_dir.subdir("share/psipred_4.01")
    temp_dir.subdir("share/psipred_4.01/data")

    # Test that it tried to install, but for some reason installation failed
    install.setup()
    out, err = capsys.readouterr()
    assert "Downloading" in out
    assert "Unpacking..." in out
    assert "Installing..." in out
    assert "\033[91mRD-MCL depends on PSIPRED, and it is not installed correctly.\033[39m\n" in out
    os.chdir(cwd)


def test_install_psipred_conda(capsys, monkeypatch):
    cwd = os.getcwd()
    temp_dir = br.TempDir()
    monkeypatch.setattr(br, "TempDir", lambda *_: temp_dir)
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)

    # Nothing is installed. User agrees to install
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    def true_or_false():
        answer_list = [False, False, False, True, True, "/Users/fakeuser/anaconda3/bin/psipred"]
        for answer in answer_list:
            yield answer

    # which returns False for the three psipred programs, true for conda,  true when checking for psipred binary, and
    # a suitable path when trying to retrieve path to psipred bin directory)
    generator = true_or_false()
    monkeypatch.setattr(shutil, "which", lambda *_, **__: next(generator))

    def generate_dirs():
        dir_list = [["bin", ["random_dir", "psipred_4.01"], ["file1", "file2"]]]
        for answer in dir_list:
            yield answer

    # Find psipred bin directory
    generator2 = generate_dirs()
    monkeypatch.setattr(os, "walk", lambda *_, **__: generator2)

    # Test installation fails when weight files are not found
    install.setup()
    out, err = capsys.readouterr()

    assert "\033[1mChecking for PSIPRED:\033[m \033[91mMissing\033[39m\n\n" in out
    assert "\033[1mCalling conda...\033[m\n" in out
    assert "\033[92mPSIPRED binary installed\033[39m\n\n" in out
    assert "\033[1mError:\033[m psi-pred data file" in out
    os.chdir(cwd)


def test_psipred_install_local_no_data(capsys, monkeypatch):
    cwd = os.getcwd()
    temp_dir = br.TempDir()
    monkeypatch.setattr(br, "TempDir", lambda *_: temp_dir)
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)

    # Nothing is installed. User wants to install.
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    # which returns False for the three psipred programs (twice) and False for two hmm programs
    monkeypatch.setattr(shutil, "which", lambda *_, **__: False)

    # Pretend to download psipred-4.01-1.tar.bz2
    #monkeypatch.setattr(urllib.request, "urlretrieve", lambda *_, **__: True)

    # Pretend psipred is installed
    temp_dir.subdir("psipred")
    bin_dir = temp_dir.subdir("psipred/bin")
    open(join(bin_dir, "psipass2"), "w").close()
    open(join(bin_dir, "psipred"), "w").close()
    open(join(bin_dir, "seq2mtx"), "w").close()

    # Test psipred found but weight files do not exist
    install.setup()
    out, err = capsys.readouterr()
    # assert False, print(out)
    assert "\033[1mChecking for PSIPRED:\033[m \033[92mFound\033[39m\n" in out
    assert "\033[1mError:\033[m psi-pred data file" in out
    os.chdir(cwd)


def test_psipred_install_local(capsys, monkeypatch):
    cwd = os.getcwd()
    temp_dir = br.TempDir()
    monkeypatch.setattr(br, "TempDir", lambda *_: temp_dir)
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)

    # Nothing is installed. User wants to install.
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    # which returns False for the three psipred programs (twice) and False for two hmm programs
    monkeypatch.setattr(shutil, "which", lambda *_, **__: False)

    # Pretend to download psipred-4.01-1.tar.bz2
    # monkeypatch.setattr(urllib.request, "urlretrieve", lambda *_, **__: True)

    # Pretend psipred is installed
    temp_dir.subdir("psipred")
    bin_dir = temp_dir.subdir("psipred/bin")
    open(join(bin_dir, "psipass2"), "w").close()
    open(join(bin_dir, "psipred"), "w").close()
    open(join(bin_dir, "seq2mtx"), "w").close()

    # Add weight files and proceed to next step
    weight_files = ["weights.dat", "weights.dat2", "weights.dat3", "weights_p2.dat",
                    "weights_s.dat", "weights_s.dat2", "weights_s.dat3"]
    data_dir = temp_dir.subdir("psipred/data")
    for file in weight_files:
        open(join(data_dir, file), "w").close()

    # Test data files found but hmmer installation fails
    install.setup()
    out, err = capsys.readouterr()
    assert "\033[1mChecking for HMMER3 programs:\033[m\n" in out
    assert "hmmbuild: \033[91mMissing\033[39m\n" in out
    assert "hmm_fwd_back: \033[91mMissing\033[39m\n" in out
    assert "\n\033[1mDownloading hmmer-3.1b2.tar.gz\033[m\n" in out
    assert "\033[1mUnpacking...\033[m\n" in out
    assert "\033[91mFailed to download HMMER3.\033" in out
    os.chdir(cwd)


def test_install_hmmer(capsys, monkeypatch):
    cwd = os.getcwd()
    temp_dir = br.TempDir()
    monkeypatch.setattr(br, "TempDir", lambda *_: temp_dir)
    monkeypatch.setattr(install, "SCRIPT_PATH", temp_dir.path)
    monkeypatch.setattr(install, "Popen", lambda *_, **__: MockPopen)

    # which does not find psipred programs, but there is a local install
    monkeypatch.setattr(shutil, "which", lambda *_, **__: False)
    temp_dir.subdir("psipred")
    bin_dir = temp_dir.subdir("psipred/bin")
    open(join(bin_dir, "psipass2"), "w").close()
    open(join(bin_dir, "psipred"), "w").close()
    open(join(bin_dir, "seq2mtx"), "w").close()

    # User wants to install everything
    monkeypatch.setattr(br, "ask", lambda *_, **__: True)

    # Pretend to download psipred-4.01-1.tar.bz2
    # monkeypatch.setattr(urllib.request, "urlretrieve", lambda *_, **__: True)

    # Add weight files and proceed to next step
    weight_files = ["weights.dat", "weights.dat2", "weights.dat3", "weights_p2.dat",
                    "weights_s.dat", "weights_s.dat2", "weights_s.dat3"]
    data_dir = temp_dir.subdir("psipred/data")
    for file in weight_files:
        open(join(data_dir, file), "w").close()

    # Create c file
    temp_dir.subdir("hmmer-3.1b2")
    src_dir = temp_dir.subdir("hmmer-3.1b2/src")
    with open(join(src_dir, "generic_fwdback.c"), 'w') as f:
        f.write("p7_gmx_Dump")

    install.setup()
    out, err = capsys.readouterr()
    assert "\n\033[1mDownloading hmmer-3.1b2.tar.gz\033[m\n" in out
    assert "\033[1mUnpacking...\033[m\n" in out
    assert "\033[1mInstalling...\033[m\n" in out
    os.chdir(cwd)
