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
from Bio.Phylo.BaseTree import TreeMixin

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
    def __init__(self, num_genes, num_taxa, seed=12345, branch_length=1.0, branch_stdev=None, drop_chance=0.0,
                 num_drops=0, duplication_chance=0.0, num_duplications=0):
        self.num_genes = num_genes
        self.num_taxa = num_taxa
        self.branch_length = branch_length

        self.rand_gen = Random(seed)  # Random seed for the hashing
        self.drop_chance = drop_chance
        self.num_drops = num_drops
        self.duplication_chance = duplication_chance
        self.num_duplications = num_duplications

        self.genes = list()  # Lists of hashes to prevent collisions
        self.taxa = list()

        for x in range(num_taxa):
            self._generate_taxa_name()
        labels = ["%s-GENE_NAME" % self.taxa[x] for x in range(num_taxa)]  # Generates a list of taxa

        self.gene_tree = Tree.randomized(num_genes, branch_length=branch_length, branch_stdev=branch_stdev)
        self.species_tree = Tree.randomized(labels, branch_length=branch_length, branch_stdev=branch_stdev)
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

        for x in range(self.num_duplications):
            leaves = tree.get_terminals()
            duplicate = self.rand_gen.random()
            if duplicate <= self.duplication_chance:
                duplicate = self.rand_gen.choice(leaves)
                duplicate.split(branch_length=self.branch_length)
                duplicate.name = ""

        for x in range(self.num_drops):
            leaves = tree.get_terminals()
            drop = self.rand_gen.random()
            if drop <= self.drop_chance:
                drop = self.rand_gen.choice(leaves)
                tree.collapse_all(drop)

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
    parser.add_argument('-nt', '--num_taxa', help='The number of taxa to be generated',
                        action='store', type=int, required=True)
    parser.add_argument('-ng', '--num_groups', help='The number of orthogroups to be generated',
                        action='store', type=int, required=True)
    parser.add_argument('-bl', '--branch_length', help='The gene tree branch length',
                        action='store', type=float, default=None)
    parser.add_argument('-bs', '--branch_stdev', help='The standard deviation of the gene tree branch length',
                        action='store', type=float, default=None)
    parser.add_argument('-drp', '--drop_chance', help='The probability of losing a species in each species tree',
                        action='store', type=float, default=0)
    parser.add_argument('-ndr', '--num_drops', help='The number of times to try dropping a species from the tree',
                        action='store', type=int, default=0)
    parser.add_argument('-dup', '--duplication_chance', help='The probability of losing a species in each species tree',
                        action='store', type=float, default=0)
    parser.add_argument('-ndp', '--num_duplications', help='The number of times to try adding a species from the tree',
                        action='store', type=int, default=0)

    in_args = parser.parse_args()
    groups = in_args.num_groups
    ntaxa = in_args.num_taxa

    in_args.branch_stdev = None if not in_args.branch_length else in_args.branch_stdev
    in_args.num_drops = 0 if in_args.drop_chance <= 0 else max(in_args.num_drops, 1)
    in_args.num_duplications = 0 if in_args.duplication_chance <= 0 else max(in_args.num_duplications, 1)

    # Tree Generation #
    generator = TreeGenerator(groups, ntaxa, branch_length=in_args.branch_length, branch_stdev=in_args.branch_stdev,
                              drop_chance=in_args.drop_chance, num_drops=in_args.num_drops,
                              duplication_chance=in_args.duplication_chance,
                              num_duplications=in_args.num_duplications)
    tree_string = str(generator)
    print(tree_string)

if __name__ == '__main__':
    main()
