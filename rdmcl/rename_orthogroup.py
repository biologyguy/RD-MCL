#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Apr 22 2015 

"""
Replace the 'group_0_XXX' designation for a given group
"""

try:
    from . import helpers as hlp

except ImportError:
    import helpers as hlp

import os
from os.path import join
import sys
import shutil
import re
from collections import OrderedDict

from buddysuite import buddy_resources as br

VERSION = hlp.VERSION
VERSION.name = "rename_orthgroup"


def prepare_clusters(ifile):
    with open(ifile, "r") as ifile:
        output = ifile.readlines()
    if output[-1] == "\n":
        del output[-1]

    for indx, line in enumerate(output):
        group_name = re.search("(^.*?)\s", line).group(1)
        line = re.sub("^.*?\s+", "", line)
        line = line.split()
        output[indx] = (group_name, line)
    output = OrderedDict(output)
    return output


def argparse_init():
    import argparse

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="rename_orthogroup", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mRename orthogroup\033[m
  What, context? Pfft...
     
  Replace default orthogroup naming with something a 
  little more meaningful to humans.
  
\033[1mUsage\033[m:
  rename_orthogroup rdmcl_dir group_name new_name [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("rdmcl_dir", help="Path to RD-MCL output directory", action="store")
    positional.add_argument("group_name", help="Group name to be changed", action="store")
    positional.add_argument("new_name", help="Replacement name (no spaces!)", action="store")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")

    parser_flags.add_argument("--force", "-f", action="store_true",
                              help="Suppress safety check and just do the rename. Be careful!!")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()
    rdmcl_dir = os.path.abspath(in_args.rdmcl_dir)
    if not os.path.isdir(rdmcl_dir):
        sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\n")
        sys.exit()

    if not os.path.isfile(join(rdmcl_dir, "final_clusters.txt")):
        sys.stderr.write("Error: The provided RD-MCL output directory "
                         "does not contain the final_clusters.txt file.\n")
        sys.exit()

    clusters = hlp.prepare_clusters(join(rdmcl_dir, "final_clusters.txt"), hierarchy=True)
    group_name = in_args.group_name
    if group_name not in clusters:
        sys.stderr.write("Error: The provided group name '%s' does not exist.\n"
                         "Here is a list of groups currently present:\n"
                         "%s\n" % (group_name, "\n".join([g for g in clusters])))
        sys.exit()

    if in_args.new_name in clusters:
        sys.stderr.write("Error: The new group name '%s' already exists, please choose a unique group name.\n" %
                         in_args.new_name)
        sys.exit()

    if " " in in_args.new_name:
        sys.stderr.write("Error: You cannot include spaces in the new name, it will mess everything else up.\n"
                         "You tried to set the new name to '%s', might I suggest '%s'?\n" %
                         (in_args.new_name, re.sub(" ", "_", in_args.new_name)))
        sys.exit()

    if "/" in in_args.new_name:
        sys.stderr.write("Error: You cannot include forward slashes (/) in the new name, it will mess up paths.\n"
                         "You tried to set the new name to '%s', might I suggest '%s'?\n" %
                         (in_args.new_name, re.sub("/", "_", in_args.new_name)))
        sys.exit()

    if not in_args.force:
        if not br.ask("You are about to permanently change all references to '{0}{2}{1}' with the new name "
                      "'{0}{3}{1}'.\n Are you sure you want to do this? y/[n]: ".format(hlp.GREEN, hlp.END,
                                                                                        group_name, in_args.new_name),
                      default=False):
            print("%sAbort!!%s" % (hlp.RED, hlp.END))
            sys.exit()

    # Modify file names
    r_files = [join("alignments", group_name), join("hmm", group_name), join("sim_scores", "%s.scores" % group_name)]
    r_files = [join(rdmcl_dir, f) for f in r_files]
    for f in r_files:
        try:
            new_path = re.sub("{0}$|({0})\.scores".format(group_name), in_args.new_name, f)
            shutil.move(f, new_path)
            print(os.path.relpath(f), "-->", os.path.relpath(new_path))
        except FileNotFoundError:
            print("Warning: '%s' does not exist." % os.path.relpath(f))

    # Modify directory names
    r_dirs = [join(rdmcl_dir, "mcmcmc", group_name)]
    for d in r_dirs:
        try:
            new_path = re.sub(r'%s($|\.scores)' % group_name, r'%s\1' % in_args.new_name, f)
            shutil.move(d, new_path)
            print(os.path.relpath(d), "-->", os.path.relpath(new_path))
        except NotADirectoryError:
            print("Warning: '%s' does not exist." % os.path.relpath(d))

    # Modify file contents
    m_files = [join(rdmcl_dir, f) for f in ["final_clusters.txt", "paralog_cliques",
                                            "manual_merge.log", "placement.log"]]
    for f in m_files:
        try:
            with open(f, "r") as ifile:
                contents = ifile.read()
            contents = re.sub(r'%s(\s|:)' % group_name, r'%s\1' % in_args.new_name, contents)
            print(os.path.relpath(f), "modified")
            with open(f, "w") as ofile:
                ofile.write(contents)
        except FileNotFoundError:
            print("Warning: '%s' does not exist." % os.path.relpath(f))


if __name__ == '__main__':
    main()
