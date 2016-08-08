#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator

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