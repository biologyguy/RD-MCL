#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 04 2016
from pyvolve import newick, evolver, model, partition

treefile = 'Ctenophores.newick'
model_tree = newick.read_tree(file=treefile)
my_model = model.Model(model_type='WAG', alpha=3, num_categories=4)
my_partition = partition.Partition(models=my_model, size=200)
evolved = evolver.Evolver(tree=model_tree, seqfile="out_seq.fa", seqfmt='nexus', partitions=my_partition)
evolved()
