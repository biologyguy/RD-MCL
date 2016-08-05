#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 05 2016

# [&R] ((A1,(B1,(C1,D1))),((A2,(B2,(C2,D2))),((A3,(B3,(C3,D3))),((A4,(B4,(C4,D4))),(A5,(B5,(C5,D5)))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),(A5,(B5,(C5,(D5,E5))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));
# [&R] ((A1,(B1,(C1,(D1,E1)))),((A2,(B2,(C2,(D2,E2)))),((A3,(B3,(C3,(D3,E3)))),((A4,(B4,(C4,(D4,E4)))),((A5,(B5,(C5,(D5,E5)))),(A6,(B6,(C6,(D6,E6)))))))));

from string import ascii_uppercase
from sys import argv
import argparse

assert len(argv) >= 3

num_taxa = int(argv[1])

groups = int(argv[2])

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

print(tree_str)
