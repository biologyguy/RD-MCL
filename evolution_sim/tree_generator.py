#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 05 2016

# [&R] ((A1,(B1,(C1,D1))),((A2,(B2,(C2,D2))),((A3,(B3,(C3,D3))),((A4,(B4,(C4,D4))),(A5,(B5,(C5,D5)))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),(A5,(B5,(C5,(D5,E5))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),(((A3,E3),(B3,(C3,D3))),(((A4,((B4,D4),(C4,E4))),(A6,(B6,(C6,(D6,E6))))),(A5,(B5,(D5,(E5,C5))))))));


import argparse
from Bio.Phylo.BaseTree import Tree
import copy
import re
import sys
import string
from string import ascii_uppercase
from random import Random
import math


def generate_perfect_tree(num_taxa, groups):
    tree_str = "[&R] "
    lefts = 0
    for group in range(groups-1):
        tree_str += "("
        lefts += 1
        for taxa in range(num_taxa-1):
            tree_str += "({0}{1},".format(ascii_uppercase[taxa], group+1)
            lefts += 1
        tree_str += "{0}{1}".format(ascii_uppercase[num_taxa-1], group+1)
        for x in range(num_taxa-1):
            tree_str += ")"
            lefts -= 1
        tree_str += ","

    for taxa in range(num_taxa-1):
        tree_str += "({0}{1},".format(ascii_uppercase[taxa], groups)
        lefts += 1
    tree_str += "{0}{1}".format(ascii_uppercase[num_taxa-1], groups)

    while lefts > 0:
        tree_str += ")"
        lefts -= 1
    tree_str += ";"
    return tree_str



class TreeGenerator:
    def __init__(self, num_genes, num_taxa, seed=12345, branch_length=1.0, branch_stdev=None):
        self.num_genes = num_genes
        self.num_taxa = num_taxa
        self.rand_gen = Random(seed)
        self.genes = list()
        self.taxa = list()
        for x in range(num_taxa):
            self._generate_taxa_name()
        labels = ["%s-GENE_NAME" % self.taxa[x] for x in range(num_taxa)]
        self.gene_tree = Tree.randomized(num_genes, branch_length=branch_length, branch_stdev=branch_stdev)
        self.species_tree = Tree.randomized(labels, branch_length=branch_length)
        self.root = self.gene_tree.clade
        self.assemble()

    def _generate_gene_name(self):
        hash_length = int(math.log(self.num_genes, len(string.ascii_uppercase)))
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
        hash_length = int(math.log(self.num_taxa, len(string.ascii_uppercase)))
        hash_length = hash_length if hash_length > 2 else 2
        while True:
            new_hash = self.rand_gen.choice(string.ascii_uppercase)
            new_hash += "".join([self.rand_gen.choice(string.ascii_lowercase) for _ in range(hash_length)])
            if new_hash in self.taxa:
                continue
            else:
                self.taxa.append(new_hash)
                return new_hash

    def _recursive_rename(self, node, gene):
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx].name = re.sub("GENE_NAME", gene, node.clades[indx].name)
            else:
                self._recursive_rename(child, gene)

    def _generate_copy(self, gene):
        tree = copy.deepcopy(self.species_tree.clade)
        self._recursive_rename(tree, gene)
        return tree

    def _recur(self, node):
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx] = self._generate_copy(self._generate_gene_name())
            else:
                self._recur(child)

    def assemble(self):
        self._recur(self.root)

    def __str__(self):
        tree_string = self.gene_tree.format("newick")
        return re.sub("\)n[0-9]*", ")", tree_string)

def main():
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

    generator = TreeGenerator(groups, ntaxa, branch_length=branch_length, branch_stdev=branch_stdev)
    tree_string = str(generator)
    print(tree_string)

if __name__ == '__main__':
    main()
