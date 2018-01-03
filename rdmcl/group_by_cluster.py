#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Aug 15 2017

"""
Output the description line for each record in each cluster, neatly organized as a group
"""

try:
    from .compare_homolog_groups import prepare_clusters
    from . import helpers
except ImportError:
    from compare_homolog_groups import prepare_clusters
    import helpers

from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import buddy_resources as br
import re
import sys
from collections import OrderedDict
import os
import argparse

VERSION = helpers.VERSION
VERSION.name = "group_by_cluster"


def make_msa(seqbuddy, aligner, trimal=()):
    """
    Create a multiple sequence alignment
    :param seqbuddy: SeqBuddy object
    :param aligner: path to alignment program
    :param trimal: List of TrimAl thresholds to try
    :return: AlignBuddy object
    """
    trimal = trimal if trimal else ["clean"]

    if len(seqbuddy) == 1:
        alignment = Alb.AlignBuddy(str(seqbuddy))
    else:
        alignment = Alb.generate_msa(Sb.make_copy(seqbuddy), aligner, quiet=True)
        ave_seq_length = Sb.ave_seq_length(seqbuddy)
        for threshold in trimal:
            align_copy = Alb.trimal(Alb.make_copy(alignment), threshold=threshold)
            cleaned_seqs = Sb.clean_seq(Sb.SeqBuddy(str(align_copy)))
            cleaned_seqs = Sb.delete_small(cleaned_seqs, 1)
            # Structured this way for unit test purposes
            if len(alignment.records()) != len(cleaned_seqs):
                continue
            elif Sb.ave_seq_length(cleaned_seqs) / ave_seq_length < 0.5:
                continue
            else:
                alignment = align_copy
                break
    return alignment


def argparse_init():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="group_by_cluster", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mGroup by Cluster\033[m
  Oooo... It's so much easier to look at now.

  Covert RD-MCL final_cluster output into sequence files,
  alignments, consensus sequences, or lists of metadata.

\033[1mUsage\033[m:
  group_by_cluster "clusters" "sequence_file" [mode] [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("clusters", action="store", help="Path to clusters file")
    positional.add_argument("sequence_file", action="store", help="Path to original sequence file")
    positional.add_argument("mode", action="store", nargs="?", default="list",
                            help="Choose the output type [list, seqs, sequences, aln, alignment, con, consensus]")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("--aligner", "-a", action="store", default="clustalo", metavar="",
                              help="Specify a multiple sequence alignment program")
    parser_flags.add_argument("--groups", "-g", action="append", nargs="+", metavar="group",
                              help="List the specific groups to process")
    parser_flags.add_argument("--max_size", "-max", action="store", type=int, metavar="",
                              help="Max cluster size to process")
    parser_flags.add_argument("--min_size", "-min", action="store", type=int, metavar="",
                              help="Min cluster size to process")
    parser_flags.add_argument("--trimal", "-trm", action="append", nargs="*", metavar="param",
                              help="Specify trimal parameters",
                              default=["gappyout", 0.5, 0.75, 0.9, 0.95, "clean"])
    parser_flags.add_argument("--strip_taxa", "-s", action="store_true",
                              help="Remove taxa prefix before searching sequence file.")
    parser_flags.add_argument("--write", "-w", action="store", metavar="", help="Specify directory to write file(s)")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()
    mode = in_args.mode.lower()
    mode = "seqs" if "sequences".startswith(mode) else mode
    mode = "aln" if "alignment".startswith(mode) else mode
    mode = "con" if "consensus".startswith(mode) else mode
    mode = "list" if "list".startswith(mode) else mode

    if mode not in ["seqs", "aln", "con", "list"]:
        Sb.br._stderr('Unrecognized mode, please select from ["seqs", "aln", "con", "list"].\n')
        sys.exit()

    if in_args.groups:
        in_args.groups = [x.lower() for x in in_args.groups[0]]
        in_args.groups = "^%s$" % "$|^".join(in_args.groups)

    cluster_file = prepare_clusters(in_args.clusters, hierarchy=True)
    seqbuddy = Sb.SeqBuddy(in_args.sequence_file)
    output = OrderedDict()

    for rank, node in cluster_file.items():
        rank = rank.split()[0]
        if in_args.groups:
            if not re.search(in_args.groups, rank):
                continue

        if in_args.min_size:
            if len(node) < in_args.min_size:
                continue

        if in_args.max_size:
            if len(node) > in_args.max_size:
                continue

        if in_args.strip_taxa:
            node = [re.sub("^.*?\-", "", x) for x in node]

        ids = "^%s$" % "$|^".join(node)
        subset = Sb.pull_recs(Sb.make_copy(seqbuddy), ids)
        subset = Sb.order_ids(subset)

        rank_output = ""
        if mode == "list":
            rank_output += rank
            for rec in subset.records:
                rec.description = re.sub("^%s" % rec.id, "", rec.description)
                rank_output += "\n%s %s" % (rec.id, rec.description)
            rank_output += "\n"

        elif mode == "seqs":
            for rec in subset.records:
                rec.description = "%s %s" % (rank, rec.description)
            rank_output += str(subset)

        elif mode in ["aln", "con"]:
            try:
                rank_output = make_msa(subset, in_args.aligner, in_args.trimal)
            except (SystemError, AttributeError) as err:
                print(err)
                sys.exit()
            rank_output.out_format = "phylip-relaxed"

        if mode == "con":
            rec = Alb.consensus_sequence(rank_output).records()[0]
            rec.id = rank
            rec.name = rank
            rec.description = ""
            rank_output.out_format = "fasta"

        output[rank] = str(rank_output)

    if not in_args.write:
        print("\n".join(data for rank, data in output.items()).strip())

    else:
        outdir = os.path.abspath(in_args.write)
        os.makedirs(outdir, exist_ok=True)
        extension = ".%s" % seqbuddy.out_format[:3] if mode == "seq" \
            else ".txt" if mode == "list" \
            else ".phy" if mode == "aln" \
            else ".fa"

        for rank, data in output.items():
            with open(os.path.join(outdir, rank + extension), "w") as ofile:
                ofile.write(data)


if __name__ == '__main__':
    main()
