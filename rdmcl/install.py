#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: install.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/RD-MCL
© license: None, this work is public domain

Description:
Ensure all RD-MCL dependencies are installed
"""

# Std library
import sys
import re
import shutil
import os
from subprocess import Popen

# My packages
from buddysuite import buddy_resources as br

# Globals
SCRIPT_PATH = os.path.abspath(os.path.dirname(__file__))


def setup():
    import urllib.request

    sys.stdout.write("\033[1mWelcome to RD-MCL!\033[m\nConfirming installation...\n\n")
    sys.stdout.write("\033[1mChecking for PSIPRED:\033[m ")
    path_install = []
    local_install = []
    programs = ["psipass2", "psipred", "seq2mtx"]
    for prog in programs:
        if shutil.which(prog):
            path_install.append(prog)
        elif os.path.isfile(os.path.join(SCRIPT_PATH, "psipred", "bin", prog)):
            local_install.append(prog)
    if sorted(list(set(path_install + local_install))) == programs:
        sys.stdout.write("\033[92mFound\033[39m\n")
        if local_install:
            path_install, local_install = False, True
        else:
            path_install, local_install = True, False
    else:
        sys.stdout.write("\033[91mMissing\033[39m\n\n")
        if br.ask("Would you like the setup script to try and install PSIPRED? [y]/n:"):
            if shutil.which("conda"):
                sys.stdout.write("\033[1mCalling conda...\033[m\n")
                path_install, local_install = True, False
                Popen("conda install -y -c biocore psipred", shell=True).wait()
            else:
                path_install, local_install = False, True
                cwd = os.getcwd()
                tmp_dir = br.TempDir()
                os.chdir(tmp_dir.path)
                if sys.platform == "darwin":
                    version = "osx-64"
                else:
                    version = "linux-64"
                url = "https://anaconda.org/biocore/psipred/4.01/download/%s/psipred-4.01-1.tar.bz2" % version
                sys.stdout.write("\n\033[1mDownloading %s binaries from %s\033[m\n" % (version, url))
                urllib.request.urlretrieve(url, "psipred-4.01-1.tar.bz2")

                sys.stdout.write("\033[1mUnpacking...\033[m\n")
                Popen("tar -xjf psipred-4.01-1.tar.bz2", shell=True).wait()

                sys.stdout.write("\033[1mInstalling...\033[m\n")
                if os.path.isdir("%s%spsipred" % (SCRIPT_PATH, os.sep)):
                    shutil.rmtree("%s%spsipred" % (SCRIPT_PATH, os.sep))
                os.makedirs("%s%spsipred" % (SCRIPT_PATH, os.sep))

                shutil.move("bin", "{0}{1}psipred{1}".format(SCRIPT_PATH, os.sep))
                shutil.move("share{0}psipred_4.01{0}data".format(os.sep),
                            "{0}{1}psipred{1}".format(SCRIPT_PATH, os.sep))
                os.chdir(cwd)

        if not shutil.which("psipred") and not os.path.isfile(os.path.join(SCRIPT_PATH, "psipred", "bin", "psipred")):
            sys.stdout.write("\033[91mRD-MCL depends on PSIPRED, and it is not installed correctly.\033[39m\n"
                             "Please see instructions at"
                             " github.com/biologyguy/RD-MCL/wiki/Installation-Guide\n\n")
            return
        else:
            sys.stdout.write("\033[92mPSIPRED binary installed\033[39m\n\n")

    # Confirm all psipred weight files are in the rdmcl directory
    weight_files = ["weights.dat", "weights.dat2", "weights.dat3", "weights_p2.dat",
                    "weights_s.dat", "weights_s.dat2", "weights_s.dat3"]
    error_msg = """\033[1mError:\033[m psi-pred data file '{0}' not found in {1}!
    Please try reinstalling PSIPRED:

       $: conda install -c biocore --no-deps --force psipred

    or build from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/

    If the problem persists, please create an issue at https://github.com/biologyguy/RD-MCL/issues
    """

    data_dir = "{0}{1}psipred{1}data".format(SCRIPT_PATH, os.sep)

    if path_install:
        os.makedirs("%s%spsipred" % (SCRIPT_PATH, os.sep), exist_ok=True)
        os.makedirs(data_dir, exist_ok=True)

        psipred_bin_dir = shutil.which("psipred").split(os.sep)[:-2]
        root, dirs, files = next(os.walk(os.sep + os.path.join(*psipred_bin_dir, "share")))
        psipred_data_dir = re.search("'(psipred.*?)'[,\]]", str(dirs)).group(1)
        psipred_data_dir = os.sep + os.path.join(*psipred_bin_dir, "share", psipred_data_dir, "data")
        for next_file in weight_files:
            if not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)) \
                    and not os.path.isfile(os.path.join(psipred_data_dir, next_file)):
                print(error_msg.format(next_file, psipred_data_dir))
                return
            elif not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)):
                shutil.copyfile(os.path.join(psipred_data_dir, next_file),
                                "{0}{1}{2}".format(data_dir, os.sep, next_file))
    elif local_install:
        for next_file in weight_files:
            if not os.path.isfile("{0}{1}{2}".format(data_dir, os.sep, next_file)):
                print(error_msg.format(next_file, data_dir))
                return

    sys.stdout.write("\033[1mChecking for HMMER3 programs:\033[m\n")
    not_installed = []
    for program in ["hmmbuild", "hmm_fwd_back"]:
        if not shutil.which(program) and not os.path.isfile(os.path.join(SCRIPT_PATH, "hmmer", program)):
            sys.stdout.write("\t\033[1m%s: \033[91mMissing\033[39m\n" % program)
            not_installed.append(program)
        else:
            sys.stdout.write("\t\033[1m%s: \033[92mFound\033[39m\n" % program)
    if not_installed:
        if br.ask("Would you like the setup script to try and install missing HMMER3 programs? [y]/n:"):
            cwd = os.getcwd()
            temp_dir = br.TempDir()
            os.chdir(temp_dir.path)
            url = "http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz"
            sys.stdout.write("\n\033[1mDownloading hmmer-3.1b2.tar.gz\033[m\n")
            urllib.request.urlretrieve(url, "hmmer-3.1b2.tar.gz")

            sys.stdout.write("\033[1mUnpacking...\033[m\n")
            Popen("tar -xzf hmmer-3.1b2.tar.gz", shell=True).wait()
            if not os.path.isdir("hmmer-3.1b2"):
                sys.stdout.write("\033[91mFailed to download HMMER3.\033[39m\nPlease see instructions at"
                                 " github.com/biologyguy/RD-MCL/wiki/Installation-Guide\n\n")
                return
            os.chdir("hmmer-3.1b2")

            sys.stdout.write("\033[1mInstalling...\033[m\n")
            Popen("./configure; make;", shell=True).wait()
            os.chdir("src")
            with open("generic_fwdback.c", "r") as ifile:
                modified = re.sub(r'(p7_gmx_Dump\(stdout, fwd, p7_DEFAULT\);)', '', ifile.read())

            with open("generic_fwdback.c", "w") as ofile:
                ofile.write(modified)

            Popen("make generic_fwdback_example", shell=True).wait()

            os.makedirs(os.path.join(SCRIPT_PATH, "hmmer"), exist_ok=True)
            try:
                if "hmm_fwd_back" in not_installed:
                    shutil.move("generic_fwdback_example", os.path.join(SCRIPT_PATH, "hmmer", "hmm_fwd_back"))
                if "hmmbuild" in not_installed:
                    shutil.move("hmmbuild", os.path.join(SCRIPT_PATH, "hmmer", "hmmbuild"))
            except FileNotFoundError:
                sys.stdout.write("\033[91mThere was a problem building the dependency HMMER3.\033[39m\n Please see"
                                 " instructions at github.com/biologyguy/RD-MCL/wiki/Installation-Guide.\n\n"
                                 "If the problem persists, please create an issue at"
                                 " https://github.com/biologyguy/RD-MCL/issues\n\n")
                return
            os.chdir(cwd)

        else:
            sys.stdout.write("\033[91mRD-MCL depends on HMMER3.\033[39m\nPlease see instructions at"
                             " github.com/biologyguy/RD-MCL/wiki/Installation-Guide\n\n")
            return

        for program in ["hmmbuild", "hmm_fwd_back"]:
            if not shutil.which(program) and not os.path.isfile(os.path.join(SCRIPT_PATH, "hmmer", program)):
                sys.stdout.write("\033[91mFailed to install HMMER3 programs.\033[39m\nPlease see instructions at"
                                 " github.com/biologyguy/RD-MCL/wiki/Installation-Guide\n\n"
                                 "If the problem persists, please create an issue at"
                                 " https://github.com/biologyguy/RD-MCL/issues\n\n")
                return
        else:
            sys.stdout.write("\033[92mHMMER3 binaries installed\033[39m\n\n")
    open(os.path.join(SCRIPT_PATH, "config.ini"), "w").close()
    sys.stdout.write("\033[1mSuccess! You're all set.\033[m\n")
    return


if __name__ == '__main__':
    setup()
