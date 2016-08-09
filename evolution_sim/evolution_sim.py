#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator
from MyFuncs import run_multicore_function

models = ['WAG'] #, 'JTT', 'LG']
taxa_range = [5] #range(2, 13)
group_range = [5] #range(2, 13)

val = .05
branch_lengths = [.1]
'''
while val <= .5:
    branch_lengths.append(val)
    val += .05
    val = round(val, 2)
'''

arguments = []
for grp in group_range:
    for tax in taxa_range:
        for mdl in models:
            for br in branch_lengths:
                arguments.append((grp, tax, mdl, br))


seed_file = 'seed_seq.raw'
with open(seed_file, 'r') as seed_io:
    seed_seq = seed_io.read()
seed_seq = seed_seq.upper().strip()


def generate(args):
    branch_stdev = 0
    alpha = None
    categories = None

    groups, taxa, model_type, branch_length = args
    out_file = "branch_outputs/{0:02d}G_{1:02d}T_{2}_{3:0.2f}B.fa".format(groups, taxa, model_type, branch_length)

    tree_file = 'branch_inputs/{0:02d}G_{1:02d}T_{2}M_{3:0.2f}B.nwk'.format(groups, taxa, model_type,
                                                                     branch_length)
    with open(tree_file, 'w') as tree_writer:
        tree_writer.write(str(TreeGenerator(taxa, groups, branch_length=branch_length,
                                            branch_stdev=branch_stdev)))
    model_tree = newick.read_tree(file=tree_file)
    my_model = model.Model(model_type=model_type, alpha=alpha, num_categories=categories)
    my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
    evolved = evolver.Evolver(tree=model_tree, partitions=my_partition)
    evolved(seqfile=out_file, seqfmt='fasta')

run_multicore_function(arguments, generate)

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