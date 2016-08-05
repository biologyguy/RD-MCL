#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition

tree_file = 'Ctenophores.newick'
seed_file = 'seed_seq.raw'
with open(seed_file, 'r') as seed_io:
    seed_seq = seed_io.read()
seed_seq = seed_seq.upper().strip()
model_tree = newick.read_tree(file=tree_file)
my_model = model.Model(model_type='WAG', alpha=3, num_categories=4)
my_partition = partition.Partition(models=my_model, root_sequence=seed_seq)
evolved = evolver.Evolver(tree=model_tree, seqfile="out_seq.fa", seqfmt='nexus', partitions=my_partition)
evolved()
