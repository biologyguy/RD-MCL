#!/usr/bin/env python3
# -*- coding: utf-8 -*-# Created on: Jan 12 2018

"""
Check how well a cluster relates to other clusters in an RD-MCL run

Only allow groups to be placed in larger groups
"""

try:
    from .compare_homolog_groups import prepare_clusters, Cluster
    from . import helpers as hlp
    from . import rdmcl
except ImportError:
    from compare_homolog_groups import prepare_clusters, Cluster
    import helpers as hlp
    import rdmcl

from buddysuite import buddy_resources as br
from buddysuite import SeqBuddy as Sb
import sys
import os
from os.path import join
from io import StringIO
import argparse
import pandas as pd
from datetime import date
from subprocess import Popen, PIPE
from multiprocessing import Lock
import scipy.stats

VERSION = hlp.VERSION
VERSION.name = "merge_orthogroups"
LOCK = Lock()


class Check(object):
    def __init__(self, rdmcl_dir):
        self.group_name = None
        self.output = None
        self.rdmcl_dir = rdmcl_dir
        self.clusters = prepare_clusters(join(self.rdmcl_dir, "final_clusters.txt"), hierarchy=True)
        self.master_clust = Cluster([seq for group, ids in self.clusters.items() for seq in ids])
        self.r_squares = pd.read_csv(join(self.rdmcl_dir, "hmm", "rsquares_matrix.csv"))
        self.fwd_scores = pd.read_csv(join(self.rdmcl_dir, "hmm", "hmm_fwd_scores.csv"))
        self.within_group_r2_df = self._prepare_within_group_r2_df()
        self.within_group_r2_dist = hlp.create_truncnorm(hlp.mean(self.within_group_r2_df.r_square),
                                                         hlp.std(self.within_group_r2_df.r_square))
        self.within_group_fwd_df = self._prepare_within_group_fwd_df()
        self.within_group_fwd_dist = scipy.stats.gaussian_kde(self.within_group_fwd_df.fwd_raw, bw_method='silverman')

        self.btw_group_r2_df = self._prepare_between_group_r2_df()
        self.btw_group_r2_dist = hlp.create_truncnorm(hlp.mean(self.btw_group_r2_df.r_square),
                                                      hlp.std(self.btw_group_r2_df.r_square))
        self.btw_group_fwd_df = self._prepare_between_group_fwd_df()
        self.btw_group_fwd_dist = scipy.stats.gaussian_kde(self.btw_group_fwd_df.fwd_raw, bw_method='silverman')

    def _prepare_within_group_r2_df(self, force=False):
        if not os.path.isfile(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv")) or force:
            sys.stderr.write("Preparing hmm/within_group_rsquares.csv...\n")
            within_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
            for g, seqs in self.clusters.items():
                if len(seqs) < 2:
                    continue
                clust_rsquares = self.r_squares.loc[(self.r_squares["rec_id1"].isin(seqs)) &
                                                    (self.r_squares["rec_id2"].isin(seqs)) &
                                                    (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()
                within_group_rsquares = within_group_rsquares.append(clust_rsquares, ignore_index=True)
            within_group_rsquares.to_csv(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv"))
        else:
            within_group_rsquares = pd.read_csv(join(self.rdmcl_dir, "hmm", "within_group_rsquares.csv"))
        return within_group_rsquares

    def _prepare_within_group_fwd_df(self, force=False):
        if not os.path.isfile(join(self.rdmcl_dir, "hmm", "within_group_fwd.csv")) or force:
            sys.stderr.write("Preparing hmm/within_group_fwd.csv...\n")
            within_group_fwd = pd.DataFrame(columns=["hmm_id", "rec_id", "fwd_raw"])
            for g, seqs in self.clusters.items():
                if len(seqs) < 2:
                    continue
                clust_fwd = self.fwd_scores.loc[(self.fwd_scores["rec_id"].isin(seqs)) &
                                                (self.fwd_scores["hmm_id"].isin(seqs))].copy()
                within_group_fwd = within_group_fwd.append(clust_fwd, ignore_index=True)
            within_group_fwd.to_csv(join(self.rdmcl_dir, "hmm", "within_group_fwd.csv"))
        else:
            within_group_fwd = pd.read_csv(join(self.rdmcl_dir, "hmm", "within_group_fwd.csv"))
        return within_group_fwd

    def _prepare_between_group_r2_df(self, force=False):
        if not os.path.isfile(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv")) or force:
            sys.stderr.write("Preparing hmm/between_group_rsquares.csv...\n")
            between_group_rsquares = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])
            i = 0
            for g1, seqs1 in self.clusters.items():
                i += 1
                if len(seqs1) < 2:
                    continue
                for g2, seqs2 in list(self.clusters.items())[i:]:
                    if len(seqs2) < 2:
                        continue
                    clust_rsquares = self.r_squares.loc[((self.r_squares["rec_id1"].isin(seqs1)) &
                                                         (self.r_squares["rec_id2"].isin(seqs2))) |
                                                        ((self.r_squares["rec_id1"].isin(seqs2)) &
                                                         (self.r_squares["rec_id2"].isin(seqs1))) &
                                                        (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()
                    between_group_rsquares = between_group_rsquares.append(clust_rsquares, ignore_index=True)
            between_group_rsquares.to_csv(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv"))
        else:
            between_group_rsquares = pd.read_csv(join(self.rdmcl_dir, "hmm", "between_group_rsquares.csv"))
        return between_group_rsquares

    def _prepare_between_group_fwd_df(self, force=False):
        file_path = join(self.rdmcl_dir, "hmm", "between_group_fwd.csv")
        if not os.path.isfile(file_path) or force:
            sys.stderr.write("Preparing hmm/between_group_fwd.csv...\n")
            between_group_fwd = pd.DataFrame(columns=["hmm_id", "rec_id", "fwd_raw"])
            i = 0
            for g1, seqs1 in self.clusters.items():
                i += 1
                if len(seqs1) < 2:
                    continue
                for g2, seqs2 in list(self.clusters.items())[i:]:
                    if len(seqs2) < 2:
                        continue
                    clust_fwd = self.fwd_scores.loc[((self.fwd_scores["hmm_id"].isin(seqs1)) &
                                                     (self.fwd_scores["rec_id"].isin(seqs2))) |
                                                    ((self.fwd_scores["hmm_id"].isin(seqs2)) &
                                                     (self.fwd_scores["rec_id"].isin(seqs1))) &
                                                    (self.fwd_scores["hmm_id"] != self.fwd_scores["rec_id"])].copy()
                    between_group_fwd = between_group_fwd.append(clust_fwd, ignore_index=True)
            between_group_fwd.to_csv(file_path)
        else:
            between_group_fwd = pd.read_csv(file_path)
        return between_group_fwd

    def check_existing_group(self, group_name):
        if group_name not in self.clusters:
            raise IndexError("Provided group name '%s' not found in named clusters: %s" %
                             (group_name, "\n".join(br.num_sorted([g for g in self.clusters]))))
        self.group_name = group_name
        query = self.clusters[group_name]
        query_score = Cluster(query, parent=self.master_clust).score()

        self.output = []
        for g, seqs in self.clusters.items():
            if g == group_name or len(seqs) < len(query):
                continue
            compare = self.r_squares.loc[((self.r_squares["rec_id1"].isin(query)) &
                                          (self.r_squares["rec_id2"].isin(seqs))) |
                                         ((self.r_squares["rec_id1"].isin(seqs)) &
                                          (self.r_squares["rec_id2"].isin(query))) &
                                         (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()

            ave, std = hlp.mean(compare.r_square), hlp.std(compare.r_square)
            upper2 = ave + (std * 2)
            upper2 = 1 if upper2 > 1 else upper2
            lower2 = ave - (std * 2)
            lower2 = 0 if lower2 < 0 else lower2
            orig_clust = Cluster(self.clusters[g], parent=self.master_clust)
            new_clust = Cluster(self.clusters[g] + query, parent=self.master_clust)
            self.output.append([g,
                                round(self.within_group_r2_dist.cdf(upper2) - self.within_group_r2_dist.cdf(lower2), 4),
                                round(self.btw_group_r2_dist.cdf(upper2) - self.btw_group_r2_dist.cdf(lower2), 4),
                                round(orig_clust.score() + query_score, 3),
                                round(new_clust.score(), 3)])
        self.output = sorted(self.output, key=lambda x: (x[1], -x[2]), reverse=True)

    @staticmethod
    def _mc_fwd_back_old_hmms(seq_chunk, args):
        try:
            hmm_scores_file, hmm_dir_path, query_file = args
            hmm_fwd_scores = pd.DataFrame(columns=["hmm_id", "rec_id", "fwd_raw"])
            for seq in seq_chunk:
                next_hmm_path = join(hmm_dir_path, "%s.hmm" % seq.id)

                fwdback_output = Popen("%s %s %s" % (rdmcl.HMM_FWD_BACK, next_hmm_path, query_file),
                                       shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].decode()
                fwd_scores_df = pd.read_csv(StringIO(fwdback_output), delim_whitespace=True,
                                            header=None, comment="#", index_col=False)
                fwd_scores_df.columns = ["rec_id", "fwd_raw", "back_raw", "fwd_bits", "back_bits"]
                fwd_scores_df["hmm_id"] = seq.id

                hmm_fwd_scores = hmm_fwd_scores.append(fwd_scores_df.loc[:, ["hmm_id", "rec_id", "fwd_raw"]],
                                                       ignore_index=True)

            hmm_fwd_scores = hmm_fwd_scores.to_csv(path_or_buf=None, header=None, index=False, index_label=False)
            with LOCK:
                with open(hmm_scores_file, "a") as ofile:
                    ofile.write(hmm_fwd_scores)
        except KeyboardInterrupt:
            pass
        return

    @staticmethod
    def _mc_run_fwd_back_new_hmm(seq_chunk, args):
        try:
            rec, out_dir, hmm_path = args
            seqbuddy = Sb.SeqBuddy(seq_chunk)
            id_hash = hlp.md5_hash("".join(sorted([seq.id for seq in seqbuddy.records])))
            seqs_file = join(out_dir, "%s.fa" % id_hash)
            seqbuddy.write(seqs_file, out_format="fasta")
            fwdback_output = Popen("%s %s %s" % (rdmcl.HMM_FWD_BACK, hmm_path, seqs_file),
                                   shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].decode()
            with LOCK:
                with open(join(out_dir, "outfile.csv"), "a") as ofile:
                    ofile.write(fwdback_output)
        except KeyboardInterrupt:
            pass
        return

    @staticmethod
    def _mc_r_squares(seq_chunks, args):
        try:
            rec, hmm_fwd_scores, output_path = args
            comparisons = pd.DataFrame(columns=["rec_id1", "rec_id2", "r_square"])

            for seq in seq_chunks:
                fwd1 = hmm_fwd_scores.loc[hmm_fwd_scores.rec_id == seq.id].sort_values(by="hmm_id")
                fwd2 = hmm_fwd_scores.loc[hmm_fwd_scores.rec_id == rec.id].sort_values(by="hmm_id")
                corr = scipy.stats.pearsonr(fwd1.fwd_raw, fwd2.fwd_raw)
                comparisons = comparisons.append(pd.DataFrame(data=[[seq.id, rec.id, corr[0]**2]],
                                                 columns=["rec_id1", "rec_id2", "r_square"]), ignore_index=True)
            comparisons = comparisons.to_csv(path_or_buf=None, header=None, index=False, index_label=False)
            with LOCK:
                with open(output_path, "a") as ofile:
                    ofile.write(comparisons)
        except KeyboardInterrupt:
            pass
        return

    def _mc_conf_inters(self, clust, args):
        try:
            g, seqs = clust
            rec, hmm_fwd_scores, output_file = args

            # Calculate R² 95% conf interval first
            compare = self.r_squares.loc[((self.r_squares["rec_id1"] == rec.id) &
                                          (self.r_squares["rec_id2"].isin(seqs))) |
                                         ((self.r_squares["rec_id1"].isin(seqs)) &
                                          (self.r_squares["rec_id2"] == rec.id)) &
                                         (self.r_squares["rec_id1"] != self.r_squares["rec_id2"])].copy()

            ave, std = hlp.mean(compare.r_square), hlp.std(compare.r_square)
            upper2 = ave + (std * 2)
            upper2 = 1.0 if upper2 > 1 else upper2
            lower2 = ave - (std * 2)
            lower2 = 0.0 if lower2 < 0 else lower2
            output = [g, len(seqs), round(lower2, 4), round(upper2, 4)]

            # Then calculate the Fwd score 95% confidence interval
            compare = hmm_fwd_scores.loc[(hmm_fwd_scores["rec_id"] == rec.id) &
                                         (hmm_fwd_scores["hmm_id"].isin(seqs))].copy()
            ave, std = hlp.mean(compare.fwd_raw), hlp.std(compare.fwd_raw)
            upper2 = ave + (std * 2)
            lower2 = ave - (std * 2)
            output += [round(lower2, 2), round(upper2, 2)]
            output = [str(i) for i in output]
            with LOCK:
                with open(output_file, "a") as ofile:
                    ofile.write("%s\n" % ",".join(output))
        except KeyboardInterrupt:
            pass
        return

    def check_new_sequence(self, rec, minimum=1):
        seqs_file = join(self.rdmcl_dir, "input_seqs.fa")
        query_file = br.TempFile()
        query_file.write(rec.format("fasta"))
        sequences = Sb.SeqBuddy(seqs_file)
        seq_chuncks = hlp.chunk_list(sequences.records, br.usable_cpu_count())
        out_dir = br.TempDir()
        hmm_path = join(self.rdmcl_dir, "hmm", "%s.hmm" % rec.id)
        br.run_multicore_function(seq_chuncks, self._mc_run_fwd_back_new_hmm, [rec, out_dir.path, hmm_path], quiet=True)

        fwd_scores_df = pd.read_csv(join(out_dir.path, "outfile.csv"), delim_whitespace=True,
                                    header=None, comment="#", index_col=False)
        fwd_scores_df.columns = ["rec_id", "fwd_raw", "back_raw", "fwd_bits", "back_bits"]
        fwd_scores_df["hmm_id"] = rec.id

        fwd_scores_df = self.fwd_scores.copy().append(fwd_scores_df.loc[:, ["hmm_id", "rec_id", "fwd_raw"]],
                                                      ignore_index=True)
        hmm_scores_file = br.TempFile()
        params = [hmm_scores_file.path, join(self.rdmcl_dir, "hmm"), query_file.path]
        seq_chuncks[-1].append(rec)
        br.run_multicore_function(seq_chuncks, self._mc_fwd_back_old_hmms, params, quiet=True)

        temp_df = pd.read_csv(hmm_scores_file.path, header=None)

        temp_df.columns = ["hmm_id", "rec_id", "fwd_raw"]
        fwd_scores_df = fwd_scores_df.append(temp_df).reset_index(drop=True)

        # Don't recalculate the entire r_squares matrix
        output_file = br.TempFile()
        br.run_multicore_function(seq_chuncks, self._mc_r_squares, [rec, fwd_scores_df, output_file.path], quiet=True)
        comparison = pd.read_csv(output_file.path, header=None)
        comparison.columns = ["rec_id1", "rec_id2", "r_square"]
        self.r_squares = self.r_squares.append(comparison, ignore_index=True)

        self.output = []
        clusters = [(g, seqs) for g, seqs in self.clusters.items() if len(seqs) >= minimum]
        output_file.clear()
        br.run_multicore_function(clusters, self._mc_conf_inters, [rec, fwd_scores_df, output_file.path], quiet=True)

        self.output = output_file.read().strip().split("\n")
        self.output = [line.strip().split(",") for line in self.output]
        self.output = sorted(self.output, key=lambda x: (float(x[4]), float(x[5])), reverse=True)
        return

    def merge(self, merge_group_name, force=False):
        merge_group = [l for l in self.output if l[0] == merge_group_name]
        if not merge_group:
            sys.stderr.write("Error: %s is not a group that %s can be merged with.\n" %
                             (merge_group_name, self.group_name))
            return
        merge_group = merge_group[0]

        do_merge = True
        if self.output[0][0] != merge_group[0] and not force:
            if not br.ask("{0}Merge Warning{1}: The group that appears to be the most\n"
                          "appropriate for {5}{2}{1} is {5}{3}{1}, but you have\n"
                          "selected {5}{4}{1}.\n"
                          "Do you wish to continue? y/[n]: ".format(hlp.RED, hlp.END, self.group_name,
                                                                    self.output[0][0], merge_group[0], hlp.GREEN),
                          default="no"):
                do_merge = False
            print()

        if merge_group[1] < 0.05 and do_merge and not force:
            if not br.ask("{0}Merge Warning{1}: Less than 5% of sequences within current\n"
                          "clusters have a similarity distribution that matches the\n"
                          "similarity distribution between {4}{2}{1} and {4}{3}{1}.\n"
                          "This makes the merge questionable.\n"
                          "Do you wish to continue? y/[n]: ".format(hlp.RED, hlp.END, self.group_name,
                                                                    merge_group_name, hlp.GREEN), default="no"):
                do_merge = False
            print()

        if merge_group[3] > merge_group[4] and do_merge and not force:
            if not br.ask("{0}Merge Warning{1}: Merging {4}{2}{1} and {4}{3}{1} will\n"
                          "reduce the combined orthogroup score, which means you will\n"
                          "increase the number of paralogs per group.\n"
                          "Do you wish to continue? y/[n]: ".format(hlp.RED, hlp.END, self.group_name,
                                                                    merge_group_name, hlp.GREEN), default="no"):
                do_merge = False
            print()

        if do_merge and \
                (force or br.ask("Last chance to abort!\n"
                                 "Merge {2}{0}{3} into {2}{1}{3}? "
                                 "y/[n]: ".format(self.group_name, merge_group_name,
                                                  hlp.GREEN, hlp.END), default="no")):

            # 1) Update clusters in self
            self.clusters[merge_group_name] = sorted(self.clusters[merge_group_name] + self.clusters[self.group_name])
            del self.clusters[self.group_name]

            # 2) Append to manual_merge.log
            with open(join(self.rdmcl_dir, "manual_merge.log"), "a") as ofile:
                ofile.write("%s %s -> %s\n" % (date.today(), self.group_name, merge_group_name))

            # 3) Rewrite final_clusters.txt
            with open(join(self.rdmcl_dir, "final_clusters.txt"), "r") as ifile:
                final_lines = ifile.readlines()
            final_lines = [l for l in final_lines if not l.startswith("%s\t" % self.group_name)]
            for indx, l in enumerate(final_lines):
                if l.startswith("%s\t" % merge_group_name):
                    final_lines[indx] = [merge_group_name, str(merge_group[4]), *self.clusters[merge_group_name]]
                    final_lines[indx] = "\t".join(final_lines[indx])
                    final_lines[indx] += "\n"
                    break
            with open(join(self.rdmcl_dir, "final_clusters.txt"), "w") as ofile:
                ofile.write("".join(final_lines))

            # 4) Update within_group_df and between_group_df files
            self._prepare_within_group_r2_df(force=True)
            self._prepare_between_group_r2_df(force=True)

            # 5) Delete group HMM
            if os.path.isfile(join(self.rdmcl_dir, "hmm", self.group_name)):
                os.remove(join(self.rdmcl_dir, "hmm", self.group_name))

            print("%sMerged!%s\n" % (hlp.GREEN, hlp.END))
        else:
            print("%sMerge aborted!%s" % (hlp.RED, hlp.END))

    def __str__(self):
        if not self.output:
            return "You must run Check.check() before printing"
        out_str = "%sTesting %s%s\n" % (hlp.BOLD, self.group_name, hlp.END)
        longest_group_name = len(sorted([g[0] for g in self.output], key=lambda x: len(x), reverse=True)[0]) + 2

        out_str += "{2}{0: <{1}}R²Within  R²Btw   OrigScore  NewScore{3}\n".format("Groups", longest_group_name,
                                                                                   hlp.UNDERLINE, hlp.END)
        for line in self.output:
            test1 = hlp.GREEN if line[1] > line[2] and line[1] >= 0.05 else hlp.RED
            test2 = hlp.GREEN if line[2] < 0.05 else hlp.RED
            test3 = hlp.GREEN if line[3] < line[4] else hlp.RED
            out_str += "{0: <{5}}{6}{1: <10}{7}{2: <8}{8}{3: <11}{4}{9}\n".format(*line, longest_group_name, test1,
                                                                                  test2, test3, hlp.DEF_FONT)
        return out_str


def argparse_init():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="merge_orthogroups", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mMerge Orthogroups\033[m
  It's not really manual curation if a computer does it for you.

  Test how well RD-MCL-generated clusters fit into larger 
  clusters from that same run.

\033[1mUsage\033[m:
  merge_orthogroups rdmcl_dir group_name [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("rdmcl_dir", action="store", help="Path to RD-MCL output directory")
    positional.add_argument("group_name", action="store", help="Name of group to test")

    # Optional commands
    parser_flags = parser.add_argument_group(title="\033[1mAvailable commands\033[m")
    parser_flags.add_argument("--merge", "-m", action="store", metavar="", help="Name of group to add to")
    parser_flags.add_argument("--force", "-f", action="store_true",
                              help="Automatically answer 'yes' to any warning messages. Use caution!")

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
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'final_clusters.txt'.\n")
        sys.exit()

    if not os.path.isfile(join(rdmcl_dir, "hmm", "rsquares_matrix.csv")):
        sys.stderr.write("Error: The provided RD-MCL output directory does not "
                         "contain the necessary file 'hmm/rsquares_matrix.csv'.\n")
        sys.exit()

    check = Check(in_args.rdmcl_dir)
    try:
        check.check_existing_group(in_args.group_name)
    except IndexError as err:
        if "Provided group name" not in str(err):
            raise err
        sys.stderr.write("%s\n" % str(err))
        sys.exit()

    print(check)

    if in_args.merge:
        check.merge(in_args.merge, in_args.force)


if __name__ == '__main__':
    main()
