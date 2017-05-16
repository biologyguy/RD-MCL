#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: May 30 2015 

"""
Convert the output from orthogroup_caller into an SVG tree of polytomies
"""
try:
    from rdmcl.compare_homolog_groups import prepare_clusters
except ImportError:
    from compare_homolog_groups import prepare_clusters


class Nexus(object):
    def __init__(self, node_list):
        self.node_list = node_list

    def print(self):
        output = "#NEXUS\n"
        output += "begin taxa;\n"
        output += "\tdimensions ntax=%s;\n" % len(self)
        output += "\ttaxlabels\n"
        for _node in self.node_list:
            for _leaf in _node.leaf_list:
                output += "\t'%s'[&!color=#%s]\n" % (_leaf.label, _leaf.color)
        output += ";\n"
        output += "end;\n\n"
        output += "begin trees;\n"
        output += "\ttree tree_1 = %s" % self.newick(self.node_list, 1)
        output += "\nend;"
        return output

    def newick(self, _nodes, depth):
        output = "("
        node_groups = {}
        for _node in _nodes:
            try:
                node_groups.setdefault(_node.rank[depth], []).append(_node)
            except IndexError:  # This happens when cliques are involved
                node_groups.setdefault(_node.rank[depth - 1], []).append(_node)

        node_groups = [_group for indx, _group in node_groups.items()]
        node_groups.sort(key=lambda x: x[0].support, reverse=True)

        for _group in node_groups:
            if len(_group) == 1:
                output += "("
                for _leaf in _group[0].leaf_list:
                    if not _leaf.support:
                        _leaf.support = 1.
                    output += "'%s':1.0," % _leaf.label
                output = "%s):1.0," % output.strip(",")

            else:
                output += self.newick(_group, depth + 1)

        output = "%s):1.0," % output.strip(",")
        if depth == 1:
            output = "%s;" % output.strip(",")
        return output

    def __len__(self):
        length = len([_leaf.label for _node in self.node_list for _leaf in _node.leaf_list])
        return length


class Node(object):
    def __init__(self, leaf_list, rank, ave_support, std_support):
        self.leaf_list = sorted(leaf_list, key=lambda x: x.support, reverse=True)
        self.rank = rank.split("_")
        self.support = ave_support
        self.std = std_support

    def __len__(self):
        return len(self.leaf_list)


class Leaf(object):
    def __init__(self, label, _support=None, _color=None):
        self.label = label
        if _support == "nan":
            _support = None
        self.support = _support
        self.color = support_color(_support)
        if _color:
            _color = [hex_string(rgb) for rgb in _color]
            red, green, blue = _color
            self.color = "%s%s%s" % (red, green, blue)


def support_color(_support):
    if not _support:
        color = "000000"

    else:
        _support = float(_support)
        if 1. < _support < 0:
            raise ValueError("Leaf support values must be between 0 and 1")

        red = hex_string(175 - (175 * _support)) if _support <= 0.5 else hex_string(0)
        green = hex_string(175 * _support) if _support >= 0.5 else hex_string(0)
        blue = hex_string(175 * (1 - (abs(0.5 - _support))))
        color = "%s%s%s" % (red, green, blue)

    return color


def hex_string(value):
    output = "0x%0.2X" % int(value)
    return output[2:]


class CtenoColors(object):
    def __init__(self):
        panx1 = ["Bab-PanxαC", "Bch-PanxαD", "Bfo-PanxαC", "Bfr-PanxαB", "BOL-PanxαE", "Cfu-PanxαA", "Cfu-PanxαB",
                 "Cfu-PanxαD", "Cfu-PanxαF", "Dgl-PanxαF", "Edu-PanxαF", "Hca-PanxαE", "Hru-PanxαD", "Lcr-PanxαC",
                 "Lla-PanxαA", "Mle-Panxα11", "Oma-PanxαD", "Pba-PanxαB", "Tin-PanxαF", "Vpa-PanxαD"]
        panx2 = ["Bab-PanxαD", "Bch-PanxαE", "Bfo-PanxαD", "Bfr-PanxαC", "BOL-PanxαC", "Cfu-PanxαC", "Dgl-PanxαG",
                 "Edu-PanxαB", "Hca-PanxαH", "Hru-PanxαE", "Lcr-PanxαJ", "Lla-PanxαB", "Mle-Panxα12", "Oma-PanxαA",
                 "Tin-PanxαE", "Vpa-PanxαG"]
        panx3 = ["Bab-PanxαB", "Bch-PanxαC", "Bfo-PanxαB", "BOL-PanxαA", "Dgl-PanxαE", "Edu-PanxαA", "Hca-PanxαB",
                 "Hru-PanxαA", "Lcr-PanxαH", "Mle-Panxα10A", "Mle-Panxα9", "Oma-PanxαC", "Tin-PanxαC", "Vpa-PanxαB"]
        panx4 = ["Bab-PanxαA", "Bch-PanxαA", "Bfo-PanxαE", "Bfr-PanxαA", "BOL-PanxαB", "Dgl-PanxαI", "Edu-PanxαE",
                 "Hca-PanxαA", "Lcr-PanxαG", "Lla-PanxαC", "Mle-Panxα3", "Pba-PanxαD", "Tin-PanxαD", "Vpa-PanxαE"]
        panx5 = ["Bab-PanxαE", "Bfo-PanxαI", "BOL-PanxαD", "Dgl-PanxαD", "Hru-PanxαB", "Lcr-PanxαD", "Mle-Panxα2",
                 "Pba-PanxαG", "Tin-PanxαB"]
        panx6 = ["Bch-PanxαB", "Bfo-PanxαJ", "BOL-PanxαF", "Hca-PanxαG", "Lcr-PanxαI", "Mle-Panxα4", "Pba-PanxαA",
                 "Vpa-PanxαA"]
        panx7 = ["Bfo-PanxαA", "Bfo-PanxαG", "Dgl-PanxαA", "Edu-PanxαD", "Hru-PanxαC", "Lcr-PanxαB", "Mle-Panxα1",
                 "Oma-PanxαB"]
        panx8 = ["Bfr-PanxαD", "Edu-PanxαH", "Hca-PanxαC", "Lcr-PanxαK", "Mle-Panxα7A", "Pba-PanxαF"]
        panx9 = ["BOL-PanxαH", "Dgl-PanxαH", "Edu-PanxαC", "Hca-PanxαF", "Mle-Panxα8", "Pba-PanxαC"]
        panx10 = ["Edu-PanxαG", "Lcr-PanxαA", "Lcr-PanxαL", "Mle-Panxα5", "Pba-PanxαE", "Vpa-PanxαF"]
        panx11 = ["Bfo-PanxαH", "Cfu-PanxαE", "Dgl-PanxαB", "Dgl-PanxαC", "Lcr-PanxαE", "Tin-PanxαA"]
        panx12 = ["Bfo-PanxαF", "Hca-PanxαD", "Mle-Panxα6"]
        panxβ = ["Hvu-PanxβA", "Hvu-PanxβB", "Hvu-PanxβC", "Hvu-PanxβD", "Hvu-PanxβE", "Hvu-PanxβF", "Hvu-PanxβG",
                 "Hvu-PanxβH", "Hvu-PanxβI", "Hvu-PanxβJ", "Hvu-PanxβK", "Hvu-PanxβL", "Hvu-PanxβM", "Hvu-PanxβN",
                 "Hvu-PanxβO"]
        orphans = ["Vpa-PanxαC", "Edu-PanxαI", "BOL-PanxαG", "Lcr-PanxαF"]

        self.genes = [panx1, panx2, panx3, panx4, panx5, panx6, panx7, panx8, panx9, panx10, panx11, panx12,
                      panxβ, orphans]

        color_map = {"panx1": (247, 116, 2), "panx2": (0, 0, 0), "panx3": (13, 73, 43), "panx4": (153, 102, 51),
                     "panx5": (103, 47, 143), "panx6": (57, 180, 74), "panx7": (8, 114, 185), "panx8": (236, 28, 36),
                     "panx9": (167, 0, 233), "panx10": (236, 30, 121), "panx11": (250, 236, 34),
                     "panx12": (20, 115, 153)}

        self.reverse_map = {}
        for indx, group in enumerate(self.genes[:-2]):
            for gene in group:
                self.reverse_map[gene] = color_map["panx%s" % (indx + 1)]

        for gene in panxβ:
            self.reverse_map[gene] = (32, 61, 255)

    def get_color(self, gene):
        try:
            return self.reverse_map[gene]
        except KeyError:
            return tuple([0, 0, 0])


def main():
    import argparse
    from buddysuite import buddy_resources as br

    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="homolog_tree_builder", formatter_class=fmt, add_help=False,
                                     usage=argparse.SUPPRESS, description='''\
\033[1mHomolog Tree Builder\033[m
  An RD-MCL output visualization tool
     
  Pass in a file containing orthogroups and it will be converted
  into an hierarchical tree!
  
\033[1mUsage\033[m:
  homolog_tree_builder "/path/to/clusters" [-options]
''')

    # Positional
    positional = parser.add_argument_group(title="\033[1mPositional argument\033[m")

    positional.add_argument("cluster_file", help="Specify a cluster file (like 'final_clusters.txt' following RD-MCL)")

    # Developer testing
    dev_flags = parser.add_argument_group(title="\033[1mDeveloper commands\033[m")
    dev_flags.add_argument("-panx", "--ctenos_panxs", action="store_true", help="Add color to Ctenophore pannexins.")

    # Misc
    misc = parser.add_argument_group(title="\033[1mMisc options\033[m")
    misc.add_argument('-h', '--help', action="help", help="Show this help message and exit")

    in_args = parser.parse_args()
    ctenos = CtenoColors()
    cluster_file = prepare_clusters(in_args.cluster_file, hierarchy=True)

    nodes = []
    for rank, node in cluster_file.items():
        leaves = []
        for leaf in node:
            if in_args.ctenos_panxs:
                leaves.append(Leaf(leaf, 1.0, ctenos.get_color(leaf)))
            else:
                leaves.append(Leaf(leaf, 1.0))

        nodes.append(Node(leaves, rank=rank, ave_support=0,
                          std_support=0))

    if in_args.ctenos_panxs:
        for node in nodes:
            temp_leaf_list = []
            for panx in ctenos.genes:
                temp_leaf_list += sorted([next_leaf for next_leaf in node.leaf_list if next_leaf.label in panx],
                                         key=lambda next_leaf: next_leaf.label)
            node.leaf_list = temp_leaf_list

    nexus = Nexus(nodes)
    print(nexus.print())


if __name__ == '__main__':
    main()