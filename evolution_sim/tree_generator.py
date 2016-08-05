#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 05 2016

# [&R] ((A1,(B1,(C1,D1))),((A2,(B2,(C2,D2))),((A3,(B3,(C3,D3))),((A4,(B4,(C4,D4))),(A5,(B5,(C5,D5)))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),(A5,(B5,(C5,(D5,E5))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),(((A3,E3),(B3,(C3,D3))),(((A4,((B4,D4),(C4,E4))),(A6,(B6,(C6,(D6,E6))))),(A5,(B5,(D5,(E5,C5))))))));

from string import ascii_uppercase
from sys import argv
import argparse
import dendropy
from dendropy.datamodel.treemodel import Tree, Node
from dendropy.datamodel.treecollectionmodel import TreeList
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.calculate import treecompare

assert len(argv) >= 3

num_taxa = int(argv[1])

groups = int(argv[2])

def generate_tree(num_taxa, groups):
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
'''
tree_str = generate_tree(num_taxa, groups)

with open("tmp.del", "w") as temp_file:
    temp_file.write(tree_str)

trees = Tree.yield_from_files(["tmp.del"], "newick")
tree = None

for _tree in trees:
    tree = _tree


print(str(tree)+";")
'''

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import copy
import re
from string import ascii_uppercase

class TreeGenerator:
    def __init__(self):
        self.letters = list(ascii_uppercase)

    def get_next_letter(self):
        return self.letters.pop()

    def recursive_rename(self, node, letter):
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx].name = re.sub("GENE_NAME", letter, node.clades[indx].name)
            else:
                self.recursive_rename(child, letter)

    def generate_copy(self, letter):
        tree = copy.deepcopy(species_tree.clade)
        self.recursive_rename(tree, letter)
        return tree

    def recur(self, node):
        for indx, child in enumerate(node):
            if child.is_terminal():
                node.clades[indx] = self.generate_copy(self.get_next_letter())
            else:
                self.recur(child)

labels = ["%sGENE_NAME" %(x+1) for x in range(num_taxa)]
gene_tree = Tree.randomized(groups)
species_tree = Tree.randomized(labels)
root = gene_tree.clade
generator = TreeGenerator()
generator.recur(root)

tree_string = gene_tree.format("newick")
tree_string = re.sub("\)n[0-9]*", ")", tree_string)
print(tree_string)