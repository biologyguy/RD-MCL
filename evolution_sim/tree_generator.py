#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 05 2016

import argparse
import copy
import math
import re
import string
from random import Random
from Bio.Phylo.BaseTree import Tree


def generate_perfect_tree(num_taxa, groups):  # Generates a perfect bipartition tree
    tree_str = "[&R] "
    lefts = 0
    for group in range(groups-1):
        tree_str += "("
        lefts += 1
        for taxa in range(num_taxa-1):
            tree_str += "({0}{1},".format(string.ascii_uppercase[taxa], group+1)
            lefts += 1
        tree_str += "{0}{1}".format(string.ascii_uppercase[num_taxa-1], group+1)
        for x in range(num_taxa-1):
            tree_str += ")"
            lefts -= 1
        tree_str += ","

    for taxa in range(num_taxa-1):
        tree_str += "({0}{1},".format(string.ascii_uppercase[taxa], groups)
        lefts += 1
    tree_str += "{0}{1}".format(string.ascii_uppercase[num_taxa-1], groups)

    while lefts > 0:
        tree_str += ")"
        lefts -= 1
    tree_str += ";"
    return tree_str


class TreeGenerator:
    def __init__(self, num_genes, num_taxa, seed=12345, branch_length=1.0, branch_stdev=None):
        self.num_genes = num_genes
        self.num_taxa = num_taxa

        self.rand_gen = Random(seed)  # Random seed for the hashing

        self.genes = list()  # Lists of hashes to prevent collisions
        self.taxa = list()

        for x in range(num_taxa):
            self._generate_taxa_name()
        labels = ["%s-GENE_NAME" % self.taxa[x] for x in range(num_taxa)]  # Generates a list of taxa

        self.gene_tree = Tree.randomized(num_genes, branch_length=branch_length, branch_stdev=branch_stdev)
        self.species_tree = Tree.randomized(labels, branch_length=branch_length)
        self.root = self.gene_tree.clade

        self._recursive_build(self.root)  # Assembles the tree

    def _generate_gene_name(self):
        hash_length = int(math.log(self.num_genes, len(string.ascii_uppercase)) + .5)
        hash_length = hash_length if hash_length > 2 else 3

        while True:
            new_hash = self.rand_gen.choice(string.ascii_uppercase)
            new_hash += "".join([self.rand_gen.choice(string.ascii_lowercase) for _ in range(hash_length)])
            if new_hash in self.genes:
                continue
            else:
                self.genes.append(new_hash)
                return new_hash

    def _generate_taxa_name(self):
        hash_length = int(math.log(self.num_taxa, len(string.ascii_uppercase)) + .5)
        hash_length = hash_length if hash_length > 2 else 2

        while True:
            new_hash = self.rand_gen.choice(string.ascii_uppercase)
            new_hash += "".join([self.rand_gen.choice(string.ascii_lowercase) for _ in range(hash_length)])
            if new_hash in self.taxa:
                continue
            else:
                self.taxa.append(new_hash)
                return new_hash

    def _recursive_rename(self, node, gene):  # Helper method for _copy_species_tree that renames nodes recursively
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx].name = re.sub("GENE_NAME", gene, node.clades[indx].name)
            else:
                self._recursive_rename(child, gene)

    def _copy_species_tree(self, gene):  # Returns a copy of the species tree with a unique gene name
        tree = copy.deepcopy(self.species_tree.clade)
        self._recursive_rename(tree, gene)
        return tree

    def _recursive_build(self, node):  # Recursively replaces the terminal nodes of the gene tree with species trees
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx] = self._copy_species_tree(self._generate_gene_name())
            else:
                self._recursive_build(child)

    def __str__(self):
        tree_string = self.gene_tree.format("newick")
        return re.sub("\)n[0-9]*", ")", tree_string)  # Removes strange node IDs that biopython includes


def main():

    # Argument Parsing #
    parser = argparse.ArgumentParser(prog="PhyloBuddy.py", usage=argparse.SUPPRESS,
                                     description='A tool for generating random ortholog trees.\nUsage:'
                                                 '\n./tree_generator.py -nt <#> -ng <#>')
    parser.add_argument('-nt', '--num_taxa', help='Specifies the number of taxa to be generated', nargs=1,
                        action='store', type=int, required=True)
    parser.add_argument('-ng', '--num_groups', help='Specifies the number of orthogroups to be generated', nargs=1,
                        action='store', type=int, required=True)
    parser.add_argument('-bl', '--branch_length', help='Specifies the gene tree branch length', nargs=1,
                        action='store', type=float, required=False)
    parser.add_argument('-bs', '--branch_stdev', help='Specifies the standard deviation of the gene tree branch length',
                        nargs=1, action='store', type=float, required=False)

    in_args = parser.parse_args()
    groups = in_args.num_groups[0]
    ntaxa = in_args.num_taxa[0]
    if in_args.branch_length:
        branch_length = in_args.branch_length[0]
    else:
        branch_length = None
    if in_args.branch_stdev and in_args.branch_length:
        branch_stdev = in_args.branch_stdev[0]
    else:
        branch_stdev = None

    # Tree Generation #
    generator = TreeGenerator(groups, ntaxa, branch_length=branch_length, branch_stdev=branch_stdev)
    tree_string = str(generator)
    print(tree_string)

if __name__ == '__main__':
    main()
