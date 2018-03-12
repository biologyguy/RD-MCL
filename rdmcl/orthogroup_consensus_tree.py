#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Aug 15 2017

"""
Create a maximum likelihood tree from consensus sequences using orthgroup sampling for support
"""

try:
    from . import helpers as hlp
    from .group_by_cluster import make_msa
except ImportError:
    import helpers as hlp
    from group_by_cluster import make_msa

from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from buddysuite import PhyloBuddy as Pb
from buddysuite import buddy_resources as br
import re
import sys
import os
from os.path import join
import argparse
import random
from multiprocessing import Lock

LOCK = Lock()
VERSION = hlp.VERSION
VERSION.name = "orthogroup_consensus_tree"


def list_features(seqbuddy):
    features = []
    for rec in seqbuddy.records:
        if not rec.features:
            raise ValueError("No features detected in %s." % rec.id)
        else:
            for feat in rec.features:
                for key, qual in feat.qualifiers.items():
                    if qual not in features:
                        features.append(qual)
    return features


def mc_bootstrap(_, args):
    all_clusters, orig_genbank, outdir = args
    seq_names = []
    recs = []
    for clust, rec_ids in all_clusters.items():
        seq_names.append(random.choice(rec_ids))
        rec = Sb.pull_recs(Sb.make_copy(orig_genbank), "^%s$" % seq_names[-1])
        Sb.rename(rec, "^%s$" % rec.records[0].id, clust)
        recs.append(rec.records[0])
    tree = Sb.SeqBuddy(recs, out_format="fasta")
    tree = Alb.generate_msa(tree, "clustalo", quiet=True)
    tree = trimal(orig_genbank, ["gappyout", 0.5, 0.75, 0.9, 0.95, "clean"], tree)
    tree = Pb.generate_tree(tree, "fasttree", quiet=True)
    clean_tree = re.sub(":[0-9]e-[0-9]+", "", str(tree))
    clean_tree = re.sub(":[0-9]+\.[0-9]+", "", clean_tree)
    clean_tree = re.sub("\)[0-9]+\.[0-9]+", ")", clean_tree)
    clean_tree = re.sub("'", "", clean_tree)
    with LOCK:
        with open(join(outdir, "bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % clean_tree)
        with open(join(outdir, "raw_bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % tree)
        with open(join(outdir, "bootstraps.txt"), "a") as ofile:
            ofile.write("%s\n" % "\t".join(sorted(seq_names)))


def trimal(seqbuddy, trimal_modes, alignment):
        # Only remove columns up to a 50% reduction in average seq length and only if all sequences are retained
        ave_seq_length = Sb.ave_seq_length(seqbuddy)
        for threshold in trimal_modes:
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

    parser = argparse.ArgumentParser(prog="orthogroup_consensus_tree", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mOrthogroup Consensus Tree Builder\033[m
  Why bootstrap on columns when you can bootstrap on whole sequences?

  Create a maximum likelihood consensus sequence tree from your orthogroups and calculate branch support by sampling
  sequences from each of those orthogroups.

\033[1mUsage\033[m:
  orthogroup_consensus_tree "rdmcl_dir" [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("rdmcl_dir", action="append", nargs="+", help="Path(s) to RD-MCL output directory(ies)")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("--aligner", "-a", action="store", default="clustalo", metavar="",
                              help="Specify a multiple sequence alignment program")
    parser_flags.add_argument("--bootstraps", "-b", action="store", type=int, default=100,
                              help="How many bootstraps do you want calculated?")
    parser_flags.add_argument("--incl_groups", "-ig", action="append", nargs="+", metavar="group",
                              help="List specific groups to process")
    parser_flags.add_argument("--excl_groups", "-eg", action="append", nargs="+", metavar="group",
                              help="List specific groups to exclude")
    # parser_flags.add_argument("--include_count", "-ic", action="store_true",
    #                          help="Add the size of each orthogroup as part of the group names")
    parser_flags.add_argument("--max_size", "-max", action="store", type=int, metavar="",
                              help="Max cluster size to process")
    parser_flags.add_argument("--min_size", "-min", action="store", type=int, metavar="",
                              help="Min cluster size to process")
    parser_flags.add_argument("--outdir", "-o", action="store", metavar="",
                              help="Specify directory to write output files to.")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def main():
    in_args = argparse_init()

    all_clusters = []
    orig_fasta = []
    orig_genbank = []
    consensus_fasta = []
    consensus_gb = []

    # Create output directory
    if in_args.outdir:
        outdir = os.path.abspath(in_args.outdir)
    else:
        outdir = join(os.getcwd(), "rdmcl_consensus_tree")
    os.makedirs(outdir, exist_ok=True)
    open(join(outdir, "bootstraps.nwk"), "w").close()
    open(join(outdir, "bootstraps.txt"), "w").close()

    # Initial processing of all sequences in the analysis (make/read consensus sequences, fetch/read prosite results)
    for rdmcl_dir in in_args.rdmcl_dir[0]:
        rdmcl_dir = os.path.abspath(rdmcl_dir)

        if not os.path.isdir(rdmcl_dir):
            sys.stderr.write("Error: The provided RD-MCL output directory does not exist.\n")
            sys.exit()

        if not os.path.isfile(join(rdmcl_dir, "final_clusters.txt")):
            sys.stderr.write("Error: The provided RD-MCL output directory '%s' does not "
                             "contain the necessary file 'final_clusters.txt'.\n" % rdmcl_dir)
            sys.exit()

        if not os.path.isfile(join(rdmcl_dir, "input_seqs.fa")):
            sys.stderr.write("Error: The provided RD-MCL output directory '%s' does not "
                             "contain the necessary file 'input_seqs.fa'.\n" % rdmcl_dir)
            sys.exit()

        orig_fasta.append(Sb.SeqBuddy(join(rdmcl_dir, "input_seqs.fa")))

        if os.path.isfile(join(rdmcl_dir, "input_seqs_psc.gb")):
            orig_genbank.append(Sb.SeqBuddy(join(rdmcl_dir, "input_seqs_psc.gb")))
            psc_records = Sb.delete_records(Sb.make_copy(orig_fasta[-1]),
                                            ["^%s$" % rec.id for rec in orig_genbank[-1].records]).records
        else:
            orig_genbank.append(None)
            psc_records = orig_fasta[-1].records

        # Get Prosite Scan results if necessary
        if psc_records:
            prosite_scan = Sb.PrositeScan(Sb.SeqBuddy(psc_records))
            prosite_scan.run()
            if orig_genbank[-1] is None:
                orig_genbank[-1] = prosite_scan.seqbuddy
            else:
                orig_genbank[-1].records += prosite_scan.seqbuddy.records
            orig_genbank[-1].write(join(rdmcl_dir, "input_seqs_psc.gb"))

        clusters = hlp.prepare_clusters(join(rdmcl_dir, "final_clusters.txt"), hierarchy=True)

        # Consensus sequences
        if os.path.isfile(join(rdmcl_dir, "consensus_seqs.fa")):
            consensus_fasta.append(Sb.SeqBuddy(join(rdmcl_dir, "consensus_seqs.fa")))
            missing_consensus = [rec.id for rec in consensus_fasta[-1].records if rec.id not in clusters]
        else:
            consensus_fasta.append(None)
            missing_consensus = list(clusters)

        for rank in missing_consensus:
            node = clusters[rank]
            ids = "^%s$" % "$|^".join(node)
            subset = Sb.pull_recs(Sb.make_copy(orig_fasta[-1]), ids)
            subset = Sb.order_ids(subset)
            try:
                rank_output = make_msa(subset, in_args.aligner)
            except (SystemError, AttributeError) as err:
                print(err)
                sys.exit()

            rec = Alb.consensus_sequence(rank_output).records()[0]
            rec.id = rank
            rec.name = rank
            rec.description = ""
            if consensus_fasta[-1] is None:
                consensus_fasta[-1] = Sb.SeqBuddy([rec], out_format="fasta")
            else:
                consensus_fasta[-1].records.append(rec)

        consensus_fasta[-1].write(join(rdmcl_dir, "consensus_seqs.fa"))

        if os.path.isfile(join(rdmcl_dir, "consensus_psc.gb")):
            consensus_gb.append(Sb.SeqBuddy(join(rdmcl_dir, "consensus_psc.gb")))
            psc_records = Sb.delete_records(Sb.make_copy(consensus_gb[-1]),
                                            ["^%s$" % rec.id for rec in consensus_fasta[-1].records]).records
        else:
            consensus_gb.append(None)
            psc_records = consensus_fasta[-1].records

        if psc_records:
            prosite_scan = Sb.PrositeScan(Sb.SeqBuddy(psc_records))
            prosite_scan.run()
            if consensus_gb[-1] is None:
                consensus_gb[-1] = prosite_scan.seqbuddy
            else:
                consensus_gb[-1].records += prosite_scan.seqbuddy.records
            consensus_gb[-1].write(join(rdmcl_dir, "consensus_psc.gb"))

        all_clusters.append(clusters)

    # Merge all of the various collections of sequences/clusters
    all_clusters = {clust: rec_ids for clusts in all_clusters for clust, rec_ids in clusts.items()}
    orig_fasta = Sb.SeqBuddy([rec for sb in orig_fasta for rec in sb.records], out_format="fasta")
    orig_genbank = Sb.SeqBuddy([rec for sb in orig_genbank for rec in sb.records], out_format="gb")
    consensus_fasta = Sb.SeqBuddy([rec for sb in consensus_fasta for rec in sb.records], out_format="fasta")
    consensus_gb = Sb.SeqBuddy([rec for sb in consensus_gb for rec in sb.records], out_format="fasta")

    # Restrict analysis to specified clusters
    group_delete = []
    recs_delete = []

    for clust, rec_ids in all_clusters.items():
        if in_args.excl_groups and clust in in_args.excl_groups[0]:
            group_delete.append(clust)
            recs_delete += rec_ids
            continue

        if in_args.min_size and len(rec_ids) < in_args.min_size:
            group_delete.append(clust)
            recs_delete += rec_ids
            continue

        if in_args.max_size and len(rec_ids) > in_args.max_size:
            group_delete.append(clust)
            recs_delete += rec_ids
            continue

        if in_args.incl_groups:
            if clust not in in_args.incl_groups[0]:
                group_delete.append(clust)
                recs_delete += rec_ids

    if group_delete:
        for clust in group_delete:
            del all_clusters[clust]
        group_delete = "^%s$" % "$|^".join(group_delete)
        Sb.delete_records(consensus_fasta, group_delete)
        Sb.delete_records(consensus_gb, group_delete)
    if recs_delete:
        recs_delete = "^%s$" % "$|^".join(recs_delete)
        Sb.delete_records(orig_fasta, recs_delete)
        Sb.delete_records(orig_genbank, recs_delete)

    # Get all features
    features = list_features(consensus_gb)
    feature_alignments = []

    # Infer consensus tree with RAxML
    cons_alignment = Alb.generate_msa(consensus_gb, "clustalo")
    Alb.trimal(cons_alignment, threshold="gappyout")
    cons_alignment.write(join(outdir, "consensus_aln.gb"))

    cons_tree = Pb.generate_tree(cons_alignment, "raxml", quiet=True)
    with open(join(outdir, "consensus_tree.nwk"), "w") as ofile:
        ofile.write(re.sub("'", "", str(cons_tree)))

    br.run_multicore_function(range(in_args.bootstraps), mc_bootstrap, func_args=[all_clusters, orig_genbank, outdir])

    support_tree = Pb.generate_tree(cons_alignment, "raxml",
                                    params="-f b -p 1234 -t %s -z %s" % (join(outdir, "consensus_tree.nwk"),
                                                                         join(outdir, "bootstraps.nwk")), quiet=True)

    support_tree.write(join(outdir, "consensus_with_support.nwk"))


if __name__ == '__main__':
    main()
