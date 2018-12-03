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
from subprocess import Popen

LOCK = Lock()
VERSION = hlp.VERSION
VERSION.name = "orthogroup_consensus_tree"


def clean_tree(tree):
    tree = re.sub(r":[0-9.]+e-[0-9]+", "", str(tree))
    tree = re.sub(r":[0-9]+\.[0-9]+", "", tree)
    tree = re.sub(r"\)[0-9]+\.[0-9]+", ")", tree)
    tree = re.sub("'", "", tree)
    return tree


def mc_hmm_bootstraps(_, args):
    all_clusters, orig_embl, outdir, feature_align, hmm_files, tree_prog, tree_params = args
    seq_names = []
    recs = []
    for clust, rec_ids in all_clusters.items():
        seq_names.append(random.choice(rec_ids))
        rec = Sb.pull_recs(Sb.make_copy(orig_embl), "^%s$" % seq_names[-1])
        try:
            Sb.rename(rec, "^%s$" % rec.records[0].id, clust)
            recs.append(rec.records[0])
        except IndexError:
            print("\nError: %s, %s" % (clust, seq_names[-1]))
            return

    seqbuddy = Sb.SeqBuddy(recs, out_format="fasta")
    temp_dir = br.TempDir()
    seqs_file = join(temp_dir.path, "seqs.fa")
    seqbuddy.write(seqs_file)
    aln_file = join(temp_dir.path, "aln.stlk")
    alignments = []
    for hmm in hmm_files:
        Popen("hmmalign --trim --amino -o %s %s %s" % (aln_file, hmm, seqs_file), shell=True).wait()
        alignments.append(Alb.AlignBuddy(aln_file).alignments[0])

    alignbuddy = Alb.AlignBuddy(alignments)
    alignbuddy = alignbuddy if len(alignbuddy) == 1 else Alb.concat_alignments(alignbuddy, group_pattern=".*")

    tree = Pb.generate_tree(alignbuddy, tree_prog, params=tree_params, quiet=True)
    with LOCK:
        with open(join(outdir, "raw_bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % tree)
        with open(join(outdir, "bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % clean_tree(tree))
        with open(join(outdir, "bootstraps.txt"), "a") as ofile:
            ofile.write("%s\n" % "\t".join(sorted(seq_names)))


def mc_bootstrap(_, args):
    all_clusters, orig_embl, outdir, feature_align, aligner, align_params, tree_prog, tree_params = args
    seq_names = []
    recs = []
    for clust, rec_ids in all_clusters.items():
        seq_names.append(random.choice(rec_ids))
        rec = Sb.pull_recs(Sb.make_copy(orig_embl), "^%s$" % seq_names[-1])
        try:
            Sb.rename(rec, "^%s$" % rec.records[0].id, clust)
            recs.append(rec.records[0])
        except IndexError:
            print("\nError: %s, %s" % (clust, seq_names[-1]))
            return

    tree = Sb.SeqBuddy(recs, out_format="fasta")
    if feature_align:
        tree = create_feature_alignment(tree, aligner, align_params)
    else:
        tree = Alb.generate_msa(tree, aligner, params=align_params, quiet=True)

    # Unless the user has specified threads, make sure RAxML doesn't try to use more resources than available
    if "raxml" in tree_prog and "-T" not in tree_params:
        tree_params += " -T 1"

    tree = Pb.generate_tree(tree, tree_prog, params=tree_params, quiet=True)
    with LOCK:
        with open(join(outdir, "raw_bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % tree)
        with open(join(outdir, "bootstraps.nwk"), "a") as ofile:
            ofile.write("%s" % clean_tree(tree))
        with open(join(outdir, "bootstraps.txt"), "a") as ofile:
            ofile.write("%s\n" % "\t".join(sorted(seq_names)))


def get_features(orig_sb, psc_file):
    recs_with_features = {}
    recs_without_features = {}

    if os.path.isfile(psc_file):
        psc_records = Sb.SeqBuddy(psc_file)
        for rec in psc_records.records:
            if rec.features:
                recs_with_features[rec.id] = rec
            else:
                recs_without_features[rec.id] = rec
        for rec in orig_sb.records:
            if rec.id not in recs_with_features and rec.id not in recs_without_features:
                recs_without_features[rec.id] = rec
    else:
        recs_without_features = {rec.id: rec for rec in orig_sb.records}

    # Get Prosite Scan results up to three times in case InterPro fails to return results for some (this can happen)
    for _ in range(3):
        if not recs_without_features:
            break

        prosite_scan = Sb.PrositeScan(Sb.SeqBuddy([rec for rec_id, rec in recs_without_features.items()]))
        prosite_scan.run()

        for rec in prosite_scan.seqbuddy.records:
            if rec.features:
                recs_with_features[rec.id] = rec
                del recs_without_features[rec.id]

    output = Sb.SeqBuddy([rec for rec_id, rec in recs_with_features.items()] +
                         [rec for rec_id, rec in recs_without_features.items()],
                         out_format="embl")

    if recs_without_features:
        sys.stderr.write("WARNING -- The following records do not have any recognized features:\n")
        for rec in recs_without_features:
            sys.stderr.write("%s\n" % rec)
        sys.stderr.write("\n")
    return output


# I think that the rest of this stuff is going to be depricated and removed
def list_features(seqbuddy):
    features = []
    for rec in seqbuddy.records:
        for feat in rec.features:
            for key, qual in feat.qualifiers.items():
                if qual[0] not in features:
                    features.append(qual[0])
    return features


def create_feature_alignment(alignment, aligner, align_params):
    features = list_features(alignment)
    feature_alignments = []

    for feature in features:
        recs_with_feat = Sb.pull_recs_with_feature(Sb.make_copy(alignment), feature)
        recs_with_feat = Sb.extract_feature_sequences(recs_with_feat, feature)
        if len(recs_with_feat) == 1:
            recs_with_feat = Alb.AlignBuddy(str(recs_with_feat), out_format="embl")
        else:
            recs_with_feat = Alb.generate_msa(recs_with_feat, aligner, params=align_params, quiet=True)
        aln_len = len(recs_with_feat.records()[0].seq)

        recs_without_feats = Sb.delete_recs_with_feature(Sb.make_copy(alignment), feature)
        for rec in recs_without_feats.records:
            rec.seq = Sb.Seq("-" * aln_len, rec.seq.alphabet)
            recs_with_feat.alignments[0].append(rec)
        feature_alignments.append(recs_with_feat.alignments[0])

    feature_alignments = Alb.AlignBuddy(feature_alignments, out_format="embl")
    Alb.concat_alignments(feature_alignments, group_pattern=".*")
    return feature_alignments


def trimal(seqbuddy, trimal_modes, alignment):
        # Only remove columns up to a 50% reduction in average seq length and only if all sequences are retained
        ave_seq_length = Sb.ave_seq_length(seqbuddy)
        for threshold in trimal_modes:
            align_copy = Alb.trimal(Alb.make_copy(alignment), threshold=threshold)
            cleaned_seqs = Sb.clean_seq(Sb.SeqBuddy(str(align_copy), out_format="embl"))
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
    parser_flags.add_argument("--align_params", "-ap", action="store", default="", metavar="",
                              help="Add any extra parameters to feed into alignment program")
    parser_flags.add_argument("--bootstraps", "-b", action="store", type=int, default=100,
                              help="How many bootstraps do you want calculated?")
    parser_flags.add_argument("--cpus", "-c", action="store", type=int, default=br.usable_cpu_count(),
                              help="Specify max CPUS")
    parser_flags.add_argument("--domain_partitions", "-dp", action="store_true",
                              help="Force alignments to PrositeScan domains")
    parser_flags.add_argument("--excl_groups", "-eg", action="append", nargs="+", metavar="group",
                              help="List specific groups to exclude")
    parser_flags.add_argument("--hmm", "-hmm", action="append", nargs="+", metavar="hmm file(s)",
                              help="Partition sequences with domain hidden Markov models")
    parser_flags.add_argument("--incl_groups", "-ig", action="append", nargs="+", metavar="group",
                              help="List specific groups to process")
    parser_flags.add_argument("--include_count", "-ic", action="store_true",
                              help="Add the size of each orthogroup as part of the group names")
    parser_flags.add_argument("--max_size", "-max", action="store", type=int, metavar="",
                              help="Max cluster size to process")
    parser_flags.add_argument("--min_size", "-min", action="store", type=int, metavar="",
                              help="Min cluster size to process")
    parser_flags.add_argument("--outdir", "-o", action="store", metavar="",
                              help="Specify directory to write output files to.")
    parser_flags.add_argument("--tree_prog", "-t", action="store", default="raxml", metavar="",
                              help="Specify a phylogenetic inference program")
    parser_flags.add_argument("--tree_prog_params", "-tp", action="store", default="", metavar="",
                              help="Add any extra parameters to feed into tree inference program")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-v', '--version', action='version', version=str(VERSION))
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    return in_args


def _prepare_rdmcl_dir(in_args, rdmcl_dir, orig_fasta, orig_embl, consensus_fasta, consensus_embl):
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
    if in_args.domain_partitions:
        orig_embl.append(get_features(orig_fasta[-1], join(rdmcl_dir, "input_seqs_psc.embl")))
    else:
        orig_embl.append(Sb.make_copy(orig_fasta[-1]))
        orig_embl[-1].out_format = "embl"

    orig_embl[-1].write(join(rdmcl_dir, "input_seqs_psc.embl"))

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

    if in_args.domain_partitions:
        consensus_embl.append(get_features(consensus_fasta[-1], join(rdmcl_dir, "consensus_psc.embl")))
    else:
        consensus_embl.append(Sb.make_copy(consensus_fasta[-1]))
        consensus_embl[-1].out_format = "embl"

    consensus_embl[-1].write(join(rdmcl_dir, "consensus_psc.embl"))
    return clusters


def main():
    in_args = argparse_init()

    all_clusters = []
    orig_fasta = []
    orig_embl = []
    consensus_fasta = []
    consensus_embl = []

    # Create output directory
    if in_args.outdir:
        outdir = os.path.abspath(in_args.outdir)
    else:
        outdir = join(os.getcwd(), "rdmcl_consensus_tree")
    os.makedirs(outdir, exist_ok=True)
    open(join(outdir, "raw_bootstraps.nwk"), "w").close()
    open(join(outdir, "bootstraps.nwk"), "w").close()
    open(join(outdir, "bootstraps.txt"), "w").close()

    # Initial processing of all sequences in the analysis (make/read consensus sequences, fetch/read prosite results)
    for rdmcl_dir in in_args.rdmcl_dir[0]:
        clusters = _prepare_rdmcl_dir(in_args, rdmcl_dir, orig_fasta, orig_embl, consensus_fasta, consensus_embl)
        all_clusters.append(clusters)

    clust_check = {}
    for clust in [c for clusts in all_clusters for c in clusts]:
        if clust in clust_check:
            sys.stderr.write("Error: Duplicate cluster detected --> '%s'.\n" % clust)
            sys.exit()
        else:
            clust_check[clust] = None

    # Merge all of the various collections of sequences/clusters
    all_clusters = {clust: rec_ids for clusts in all_clusters for clust, rec_ids in clusts.items()}
    orig_fasta = Sb.SeqBuddy([rec for sb in orig_fasta for rec in sb.records], out_format="fasta")
    orig_embl = Sb.SeqBuddy([rec for sb in orig_embl for rec in sb.records], out_format="embl")
    consensus_fasta = Sb.SeqBuddy([rec for sb in consensus_fasta for rec in sb.records], out_format="fasta")
    consensus_embl = Sb.SeqBuddy([rec for sb in consensus_embl for rec in sb.records], out_format="embl")

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
        Sb.delete_records(consensus_embl, group_delete)

    if recs_delete:
        recs_delete = "^%s$" % "$|^".join(recs_delete)
        Sb.delete_records(orig_fasta, recs_delete)
        Sb.delete_records(orig_embl, recs_delete)

    # Hash cluster names, to prevent any weirdness
    hash_dict = {}
    for clust in all_clusters:
        clust_hash = hlp.md5_hash(clust)
        hash_dict[clust_hash] = clust
        consensus_embl = Sb.rename(consensus_embl, "^%s$" % clust, clust_hash)

    for clust_hash, clust in hash_dict.items():
        all_clusters[clust_hash] = all_clusters[clust]
        del all_clusters[clust]

    if in_args.include_count:
        for clust, rec_ids in all_clusters.items():
            hash_dict[clust] += "--%s" % len(rec_ids)

    # Infer consensus
    if in_args.domain_partitions:
        cons_alignment = create_feature_alignment(consensus_embl, in_args.aligner, in_args.align_params)
    elif in_args.hmm:
        hmm_files = in_args.hmm[0]
        for hmm in hmm_files:
            if not os.path.isfile(hmm):
                sys.stderr.write("Error: The provided HMM path '%s' is not a file.\n" % hmm)
                sys.exit()
        in_args.hmm = hmm_files
        in_args.aligner = "hmmalign"
        in_args.align_params = " --trim --amino"
        temp_dir = br.TempDir()
        seqs_file = join(temp_dir.path, "seqs.fa")
        consensus_embl.write(seqs_file, out_format="fasta")
        aln_file = join(temp_dir.path, "aln.stlk")
        alignments = []
        for hmm in hmm_files:
            Popen("hmmalign --trim --amino -o %s %s %s" % (aln_file, hmm, seqs_file), shell=True).wait()
            alignments.append(Alb.AlignBuddy(aln_file).alignments[0])
        cons_alignment = Alb.AlignBuddy(alignments)
        cons_alignment = cons_alignment if len(cons_alignment) == 1 else Alb.concat_alignments(cons_alignment,
                                                                                               group_pattern=".*")
    else:
        cons_alignment = Alb.generate_msa(consensus_embl, "clustalo", params=in_args.align_params, quiet=True)

    cons_alignment.write(join(outdir, "consensus_aln.embl"), out_format="embl")

    sys.stderr.write("""
Multiple sequence Alignement:
    $: %s %s

Phylogenetic inference (%s bootstraps):
    $: %s %s

""" % (in_args.aligner, in_args.align_params, in_args.bootstraps, in_args.tree_prog, in_args.tree_prog_params))

    cons_tree = Pb.generate_tree(cons_alignment, in_args.tree_prog, params=in_args.tree_prog_params, quiet=True)
    with open(join(outdir, "consensus_tree.nwk"), "w") as ofile:
        ofile.write(re.sub("'", "", str(cons_tree)))

    if in_args.hmm:
        func_args = [all_clusters, orig_embl, outdir, in_args.domain_partitions, in_args.hmm,
                     in_args.tree_prog, in_args.tree_prog_params]

        br.run_multicore_function(range(in_args.bootstraps), mc_hmm_bootstraps,
                                  func_args=func_args, max_processes=in_args.cpus)
    else:
        func_args = [all_clusters, orig_embl, outdir, in_args.domain_partitions,
                     in_args.aligner, in_args.align_params, in_args.tree_prog, in_args.tree_prog_params]

        br.run_multicore_function(range(in_args.bootstraps), mc_bootstrap,
                                  func_args=func_args, max_processes=in_args.cpus)

    support_tree = Pb.generate_tree(cons_alignment, "raxml",
                                    params="-f b -p 1234 -t %s -z %s" % (join(outdir, "consensus_tree.nwk"),
                                                                         join(outdir, "bootstraps.nwk")), quiet=True)

    support_tree.write(join(outdir, "consensus_with_support.nwk"))

    for next_file in ["consensus_tree.nwk", "bootstraps.nwk", "consensus_with_support.nwk",
                      "consensus_aln.embl", "raw_bootstraps.nwk"]:
        for clust_hash, clust_name in hash_dict.items():
            with open(join(outdir, next_file), "r") as ifile:
                data = re.sub(r'%s' % clust_hash, r'%s' % clust_name, ifile.read())

            data = re.sub("'", "", data) if next_file.endswith("nwk") else data
            with open(join(outdir, next_file), "w") as ofile:
                ofile.write(data)
    return


if __name__ == '__main__':
    main()
