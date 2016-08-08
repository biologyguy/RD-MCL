#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition
from tree_generator import TreeGenerator

tree_file = 'tmp.nwk' #'Ctenophores.newick'
with open(tree_file, 'w') as tree_writer:
    tree_writer.write(str(TreeGenerator(10, 10, branch_length=.1)))

seed_file = 'seed_seq.raw'
with open(seed_file, 'r') as seed_io:
    seed_seq = seed_io.read()
seed_seq = seed_seq.upper().strip()
model_tree = newick.read_tree(file=tree_file)
my_model = model.Model(model_type='WAG', alpha=3, num_categories=4)
my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
evolved = evolver.Evolver(tree=model_tree, seqfile="out_seq.fa", seqfmt='nexus', partitions=my_partition)
evolved()
