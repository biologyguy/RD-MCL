#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator
from MyFuncs import run_multicore_function
import os

group_range = [5]  # range(2, 13)
taxa_range = [5]  # range(2, 13)
models = ['WAG']  # , 'JTT', 'LG']

val = .05
branch_lengths = [.1]
'''
while val <= .5:
    branch_lengths.append(val)
    val += .05
    val = round(val, 2)
'''

branch_stdevs = [0]
alphas = [None]
category_range = [None]
drop_chances = [0]
num_drops_range = [0]
duplication_chances = [0]
num_duplications_range = [0]

directory = 'DIRECTORY NAME HERE'  # Please name your run
os.makedirs('{0}_inputs'.format(directory), exist_ok=True)
os.makedirs('{0}_outputs'.format(directory), exist_ok=True)
with open("{0}_runs.txt".format(directory), "w") as runs_file:
    runs_file.write("runid\tgroup\ttaxa\tmodel\tbranch\tstdv\talpha\tcat#\tdrop\tdrop#\tdup\t"
                    "dup#\n")

seed_file = 'seed_seq.raw'
with open(seed_file, 'r') as seed_io:
    seed_seq = seed_io.read()
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


def generate(args, dirname):
    dirname = dirname[0]
    groups, taxa, model_type, branch_length, branch_stdev, alpha, categories, drop_chance, num_drops,\
        duplication_chance, num_duplications, run_id = args
    out_file = "{0}_outputs/run{1:04d}.fa".format(dirname, run_id)

    tree_file = "{0}_inputs/run{1:04d}.fa".format(dirname, run_id)

    with open("{0}_inputs/runs.txt".format(dirname), "a") as runs_file:
        runs_file.write("run{0:04d}\t{1:02d}\t{2:02d}\t{3}\t{4:.2f}\t{5:.2f}\t{6}\t{7}\t"
                        "{8:.2f}\t{9:02d}\t{10:.2f}\t{11:02d}\n".format(run_id, groups, taxa, model_type, branch_length,
                                                                        branch_stdev, alpha, categories, drop_chance,
                                                                        num_drops, duplication_chance,
                                                                        num_duplications))

    with open(tree_file, 'w') as tree_writer:
        tree_writer.write(str(TreeGenerator(taxa, groups, branch_length=branch_length, branch_stdev=branch_stdev,
                                            drop_chance=drop_chance, num_drops=num_drops,
                                            duplication_chance=duplication_chance, num_duplications=num_duplications)))
    model_tree = newick.read_tree(file=tree_file)
    my_model = model.Model(model_type=model_type, alpha=alpha, num_categories=categories)
    my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
    evolved = evolver.Evolver(tree=model_tree, partitions=my_partition)
    evolved(seqfile=out_file, seqfmt='fasta')

run_multicore_function(arguments, generate, func_args=[directory])

'''
groups = 3
taxa = 4

alpha = 3
categories = 4
branch_length = .05
branch_stdev = .1
model_type = 'WAG'

out_file = "out_seq.fa"

tree_file = 'tmp.nwk'  # 'Ctenophores.newick'
with open(tree_file, 'w') as tree_writer:
    tree_writer.write(str(TreeGenerator(taxa, groups, branch_length=branch_length, branch_stdev=branch_stdev)))

seed_file = 'seed_seq.raw'
with open(seed_file, 'r') as seed_io:
    seed_seq = seed_io.read()
seed_seq = seed_seq.upper().strip()
model_tree = newick.read_tree(file=tree_file)
my_model = model.Model(model_type=model_type, alpha=alpha, num_categories=categories)
my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
evolved = evolver.Evolver(tree=model_tree, partitions=my_partition)
evolved(seqfile=out_file, seqfmt='fasta')
'''