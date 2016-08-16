#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator
from multiprocessing import Process, SimpleQueue, Pipe
from MyFuncs import run_multicore_function, TempFile, DynamicPrint
import argparse
import subprocess
import json
import os
import re
import sqlite3
from buddysuite import SeqBuddy as Sb


def broker_func(queue):
    sqlite_file = "{0}.sqlite".format(run_name)
    if os.path.exists(sqlite_file):
        os.remove(sqlite_file)
    connection = sqlite3.connect(sqlite_file)
    cursor = connection.cursor()
    cursor.execute("CREATE TABLE data_table (runid INTEGER PRIMARY KEY)")
    cols = [("group", "INTEGER"), ("taxa", "INTEGER"), ("model", "TEXT"), ("branch", "REAL"), ("stdev", "TEXT"),
            ("alpha", "REAL"), ("catnum", "INTEGER"), ("drop", "REAL"), ("dropnum", "INTEGER"), ("dup", "REAL"),
            ("dupnum", "INTEGER"), ("tree", "TEXT"), ("groups_file", "TEXT"), ("alignment", "TEXT")]
    for pair in cols:
        cursor.execute("ALTER TABLE data_table ADD COLUMN '{0}' {1}".format(pair[0], pair[1]))

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
    groups, taxa, model_type, branch_length, branch_stdev, alpha, categories, drop_chance, num_drops,\
        duplication_chance, num_duplications, run_id = args
    generator = TreeGenerator(taxa, groups, branch_length=branch_length, branch_stdev=branch_stdev,
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

    fields = [("group", groups), ("taxa", taxa), ("model", model_type), ("branch", branch_length),
              ("stdev", branch_stdev), ("alpha", alpha), ("catnum", categories), ("drop", drop_chance),
              ("dropnum", num_drops), ("dup", duplication_chance), ("dupnum", num_duplications),
              ("tree", tree_string), ("groups_file", groups_string), ("alignment", alignment)]

    push("INSERT INTO data_table (runid) VALUES ({})".format(run_id))
    for field in fields:
        command = "UPDATE data_table SET '{0}'=('{1}') WHERE runid=({2})".format(field[0], field[1], run_id)
        push(command)

    return


def rd_mcl():
    print()
    printer = DynamicPrint()
    for run in range(num_runs):
        tmp_file = TempFile()
        msa = fetch("SELECT (alignment) FROM data_table WHERE runid='{0}'".format(run))
        tmp_file.write(msa)
        outdir = "{0}_{1:05d}".format(run_name, run)
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


def parse_input(input_str, datatype):
    if input_str is None or input_str is [None]:
        return [None]
    input_str = re.split(' ', input_str)
    if len(input_str) == 1:
        return [datatype(input_str[0])]
    elif len(input_str) == 2:
        return make_range(datatype(input_str[0]), datatype(input_str[1]), 1)
    elif len(input_str) == 3:
        return make_range(datatype(input_str[0]), datatype(input_str[1]), datatype(input_str[2]))
    else:
        raise AttributeError("Too many arguments provided: %s" % input_str)


if __name__ == '__main__':
    # Argument Parsing #
    parser = argparse.ArgumentParser(prog="evolution_sim.py", usage=argparse.SUPPRESS,
                                     description='A tool for generating orthogroups and running RD-MCL on them.\nUsage:'
                                                 '\n./evolution_sim.py seed_file run_name -nt <#> -ng <#> ...')
    parser.add_argument('seed_file', metavar='seqfile', type=str, help="The seed sequence file")
    parser.add_argument('run_name', metavar='name', type=str, help="A name for your run")
    parser.add_argument('-nt', '--num_taxa',
                        help='The number of taxa to be generated. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]")',
                        action='store', required=True)
    parser.add_argument('-ng', '--num_groups',
                        help='The number of orthogroups to be generated. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', required=True)
    parser.add_argument('-bl', '--branch_length',
                        help='The gene tree branch length. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default=None)
    parser.add_argument('-bs', '--branch_stdev',
                        help='The standard deviation of the gene tree branch length. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default=None)
    parser.add_argument('-drp', '--drop_chance',
                        help='The probability of losing a species in each species tree. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default='0')
    parser.add_argument('-ndr', '--num_drops',
                        help='The number of times to try dropping a species from the tree. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default='0')
    parser.add_argument('-dup', '--duplication_chance',
                        help='The probability of losing a species in each species tree. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default='0')
    parser.add_argument('-ndp', '--num_duplications',
                        help='The number of times to try adding a species from the tree. '
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default='0')
    parser.add_argument('-mdl', '--models',
                        help='The number of times to try adding a species from the tree. '
                             'args: list of model names "mdl1 mdl2 mdl3". Options: WAG, JTT, LG',
                        action='store', default='WAG')
    parser.add_argument('-alp', '--alpha',
                        help='The shape parameter of the gamma distribution.'
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default=None)
    parser.add_argument('-cat', '--categories',
                        help='The number of gamma categories.'
                             '1-3 args: "value" or "lower_bound upper_bound step[default=1]"',
                        action='store', default=None)
    parser.add_argument('-skp', '--skip_rdmcl',
                        help='Don\'t run RD-MCL on the generated data.', action='store_true', default=False)
    parser.add_argument('-arg', '--rdmcl_args',
                        help='Arguments to be passed to RD-MCL.',
                        action='store', default='')

    in_args = parser.parse_args()

    call_rdmcl = not in_args.skip_rdmcl
    run_name = in_args.run_name  # Please name your run
    rdmcl_flags = in_args.rdmcl_args

    group_range = parse_input(in_args.num_groups, int)
    taxa_range = parse_input(in_args.num_taxa, int)
    models = in_args.models.split(' ')

    branch_lengths = parse_input(in_args.branch_length, float)

    branch_stdevs = parse_input(in_args.branch_stdev, float)
    alphas = parse_input(in_args.alpha, float)
    category_range = parse_input(in_args.categories, int)
    drop_chances = parse_input(in_args.drop_chance, float)
    num_drops_range = parse_input(in_args.num_drops, int)
    duplication_chances = parse_input(in_args.duplication_chance, float)
    num_duplications_range = parse_input(in_args.num_duplications, int)

    seed_file = in_args.seed_file
    assert os.path.exists(seed_file)

    with open(seed_file, 'r') as seed_io:
        seed_seq = seed_io.read()
    seed_seq = str(Sb.clean_seq(Sb.SeqBuddy(seed_seq, out_format='raw')))
    seed_seq = seed_seq.upper().strip()

    # ugly-ass loop
    arguments = []
    runid = 0
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
                                                arguments.append((grp, tax, mdl, br, stdv, alp, cat, drp, ndr, dup, ndp,
                                                                  runid))
                                                runid += 1
    num_runs = runid

    broker_queue = SimpleQueue()
    broker = Process(target=broker_func, args=[broker_queue])
    broker.start()

    run_multicore_function(arguments, generate)

    if call_rdmcl:
        rd_mcl()

    while not broker_queue.empty():
        pass

    broker.terminate()
