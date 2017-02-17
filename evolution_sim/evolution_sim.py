#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator
from multiprocessing import Process, SimpleQueue, Pipe
from MyFuncs import run_multicore_function, TempFile, DynamicPrint, TempDir
import argparse
import shutil
import subprocess
import json
import os
import re
import sys
import sqlite3
from buddysuite import SeqBuddy as Sb
import numpy as np


def broker_func(queue):
    sqlite_file = "{0}/{0}.sqlite".format(output_dir)
    connection = sqlite3.connect(sqlite_file)
    cursor = connection.cursor()
    try:
        cols = [("num_groups", "INTEGER"), ("num_taxa", "INTEGER"), ("model", "TEXT"), ("branch", "REAL"), ("stdev", "TEXT"),
                ("alpha", "REAL"), ("catnum", "INTEGER"), ("drop_chance", "REAL"), ("dropnum", "INTEGER"), ("dup", "REAL"),
                ("dupnum", "INTEGER"), ("tree", "TEXT"), ("groups_file", "TEXT"), ("alignment", "TEXT")]

        cursor.execute("CREATE TABLE data_table (%s)" % ", ".join(["%s %s" % (a, b) for a, b in cols]))
    except sqlite3.OperationalError as err:
        if "table data_table already exists" not in str(err):
            raise err

    while True:
        if not queue.empty():
            output = queue.get()
            if output[0] == 'push':
                command = output[1]
                cursor.execute(command)
                connection.commit()
            elif output[0] == 'fetch':
                command, pipe = output[1:3]
                cursor.execute(command)
                response = cursor.fetchone()
                if response is None or len(response) == 0:
                    response = None
                else:
                    response = response[0]
                pipe.send(json.dumps(response))
            else:
                raise RuntimeError("Invalid instruction for broker. Must be either put or fetch, not %s." % output[0])


def push(command):
    broker_queue.put(('push', command))


def fetch(command):
    recvpipe, sendpipe = Pipe(False)
    broker_queue.put(('fetch', command, sendpipe))
    response = json.loads(recvpipe.recv())
    return response


def generate(args):
    num_groups, num_taxa, model_type, branch_length, branch_stdev, alpha, categories, drop_chance, num_drops,\
        duplication_chance, num_duplications = args
    generator = TreeGenerator(num_taxa, num_groups, branch_length=branch_length, branch_stdev=branch_stdev,
                              drop_chance=drop_chance, num_drops=num_drops, duplication_chance=duplication_chance,
                              num_duplications=num_duplications)

    tree_string = str(generator)
    groups_string = []
    for group in generator.groups():
        for seq in group:
            groups_string.append(seq + '\t')
        groups_string.append('\n')
    groups_string = "".join(groups_string)

    model_tree = newick.read_tree(tree=tree_string)
    my_model = model.Model(model_type=model_type, alpha=alpha, num_categories=categories)
    my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
    evolved = evolver.Evolver(tree=model_tree, partitions=my_partition)

    alignment_file = TempFile()

    evolved(seqfmt='fasta', seqfile=alignment_file.path)
    with open(alignment_file.path, 'r') as aln_data:
        alignment = aln_data.read()

    fields = ["num_groups", "num_taxa", "model", "branch", "stdev", "alpha", "catnum", "drop_chance", "dropnum", "dup", "dupnum",
              "tree", "groups_file", "alignment"]

    values = [num_groups, num_taxa, "'%s'" % model_type, branch_length, "'%s'" % branch_stdev, alpha, categories, drop_chance, num_drops,
              duplication_chance, num_duplications, "'%s'" % tree_string, "'%s'" % groups_string, "'%s'" % alignment]
    values = [str(value) for value in values]

    push("INSERT INTO data_table (%s) VALUES (%s)" % (",".join(fields), ", ".join(values)))
    return


def rd_mcl():  # This is broken, because num_runs...
    print()
    printer = DynamicPrint()
    for run in range(num_runs):
        tmp_file = TempFile()
        msa = fetch("SELECT (alignment) FROM data_table WHERE runid='{0}'".format(run))
        tmp_file.write(msa)
        outdir = "{0}/rdmcl_{1:05d}".format(output_dir, run)
        printer.write("Processing run {0} of {1}...".format(run, num_runs-1))
        subprocess.run("../rdmcl/rdmcl.py {0} {1} {2}".format(tmp_file.path, outdir, rdmcl_flags), shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print()


def make_range(lower, upper, step):
    range_list = []
    val = lower
    while val < upper:
        range_list.append(val)
        val += step
    return range_list


def make_range_from_inargs(_input):
    if not _input:
        return [None]
    if len(_input) in [1, 2]:
        return _input
    elif len(_input) == 3:
        _range = np.arange(_input[0], _input[1] + 0.0000001, _input[2])  # The added value is to force 'inclusive'
        return list(_range)
    else:
        raise AttributeError("Too many arguments provided: %s" % _input)


if __name__ == '__main__':
    # Argument Parsing #
    parser = argparse.ArgumentParser(prog="evolution_sim.py", usage=argparse.SUPPRESS,
                                     description='A tool for generating orthogroups and running RD-MCL on them.\nUsage:'
                                                 '\n./evolution_sim.py seed_file output_dir -nt <#> -ng <#> ...')
    parser.add_argument('seed_file', metavar='seqfile', type=str, help="The seed sequence file (most formats are fine)")
    parser.add_argument('output_dir', metavar='name', type=str, help="A name for your run")

    parser.add_argument('-nt', '--num_taxa', nargs='+', type=int, required=True,
                        help='The number of taxa to be generated. Specify number or range')
    parser.add_argument('-ng', '--num_groups', nargs='+', type=int, required=True,
                        help='The number of orthogroups to be generated. Specify number or range')

    parser.add_argument('-bl', '--branch_length', nargs="+", type=float, default=[1],
                        help='The average gene tree branch length as substitutions per position')
    parser.add_argument('-bs', '--branch_stdev', nargs="+", type=float, default=[0.2],
                        help='The standard deviation around the gene tree branch length')
    parser.add_argument('-drp', '--drop_chance', nargs='+', type=float, default=[0],
                        help='The probability of losing a species in each species tree')
    parser.add_argument('-ndr', '--num_drops', nargs='+', type=int, default=[0],
                        help='The number of times to try dropping a species from the tree')
    parser.add_argument('-dup', '--duplication_chance', nargs='+', type=float, default=[0],
                        help='The probability of losing a species in each species tree')
    parser.add_argument('-ndp', '--num_duplications', nargs='+', type=int, default=[0],
                        help='The number of times to try adding a species from the tree')
    parser.add_argument('-mdl', '--models', nargs='+', default=['WAG'], choices=['WAG', 'JTT', 'LG'],
                        help='What evolutionary model do you want to apply?')
    parser.add_argument('-alp', '--alpha', nargs='+', type=float, default=[2],
                        help='The shape parameter of the gamma distribution')
    parser.add_argument('-cat', '--categories', nargs='+', type=int, default=[3],
                        help='The number of gamma categories.')

    parser.add_argument('-skp', '--skip_rdmcl',
                        help='Don\'t run RD-MCL on the generated data.', action='store_true', default=False)
    parser.add_argument('-arg', '--rdmcl_args',
                        help='Arguments to be passed to RD-MCL.',
                        action='store', default='')

    in_args = parser.parse_args()

    call_rdmcl = not in_args.skip_rdmcl
    output_dir = in_args.output_dir  # Please name your run
    rdmcl_flags = in_args.rdmcl_args

    taxa_range = sorted(in_args.num_taxa)
    taxa_range = taxa_range if len(taxa_range) == 1 else list(range(taxa_range[0], taxa_range[1] + 1))

    group_range = sorted(in_args.num_groups)
    group_range = group_range if len(group_range) == 1 else list(range(group_range[0], group_range[1] + 1))

    branch_lengths = make_range_from_inargs(in_args.branch_length)
    branch_stdevs = make_range_from_inargs(in_args.branch_stdev)

    drop_chances = make_range_from_inargs(in_args.drop_chance)
    if drop_chances[-1] >= 1:
        raise ValueError("A drop_chance parameter was >= 1.0")

    num_drops_range = sorted(in_args.num_drops)
    num_drops_range = num_drops_range if len(num_drops_range) == 1 \
        else list(range(num_drops_range[0], num_drops_range[1] + 1))

    duplication_chances = make_range_from_inargs(in_args.duplication_chance)

    num_duplications_range = sorted(in_args.num_duplications)
    num_duplications_range = num_duplications_range if len(num_duplications_range) == 1 \
        else list(range(num_duplications_range[0], num_duplications_range[1] + 1))

    models = in_args.models
    alphas = make_range_from_inargs(in_args.alpha)

    category_range = sorted(in_args.categories)
    category_range = category_range if len(category_range) == 1 \
        else list(range(num_drops_range[0], num_drops_range[1] + 1))

    seed_file = in_args.seed_file
    assert os.path.exists(seed_file)

    with open(seed_file, 'r') as seed_io:
        seed_seq = seed_io.read()
    seed_seq = str(Sb.clean_seq(Sb.SeqBuddy(seed_seq, out_format='raw')))
    seed_seq = seed_seq.upper().strip()

    # ugly-ass loop
    arguments = []
    for grp in group_range:
        for tax in taxa_range:
            for mdl in models:
                for br in branch_lengths:
                    for stdv in branch_stdevs:
                        for alp in alphas:
                            for cat in category_range:
                                for drp in drop_chances:
                                    for ndr in num_drops_range:
                                        for dup in duplication_chances:
                                            for ndp in num_duplications_range:
                                                arguments.append((grp, tax, mdl, br, stdv, alp, 
                                                                  cat, drp, ndr, dup, ndp))

    broker_queue = SimpleQueue()
    broker = Process(target=broker_func, args=[broker_queue])
    broker.daemon = True
    broker.start()

    os.makedirs(output_dir, exist_ok=True)

    run_multicore_function(arguments, generate)

    with open("site_rates_info.txt", "a") as ofile:
        ofile.write("\n")
    os.remove("site_rates_info.txt")
    os.remove("site_rates.txt")

    # Broken
    #if call_rdmcl:
    #    rd_mcl()

    while not broker_queue.empty():
        pass

    broker.terminate()
