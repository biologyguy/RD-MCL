#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 1 2017 

"""
DESCRIPTION OF PROGRAM
"""

import sys
from collections import OrderedDict
from buddysuite import SeqBuddy as Sb
import re


def add_tax(rec_id, taxonomy, tax_dict):
    tax_dict.setdefault("rec_ids", [])
    tax_dict["rec_ids"].append(rec_id)
    if len(taxonomy) >= 1:
        tax_dict.setdefault(taxonomy[0], OrderedDict())
        tax_dict[taxonomy[0]] = add_tax(rec_id, taxonomy[1:], tax_dict[taxonomy[0]])
    return tax_dict


def to_string(tax_dict, level=0, max_level=5, print_ids=False):
    output = "Total records: %s\n\n" % len(tax_dict["rec_ids"]) if level == 0 else ""
    if not max_level or level <= max_level:
        del tax_dict["rec_ids"]
        tax_dict = sorted([(taxon, sub_dict) for taxon, sub_dict in tax_dict.items()],
                          key=lambda x: len(x[1]))
        for taxon, sub_dict in tax_dict:
            output += "{0}{1}    {2}\n".format(" |" * level, taxon, len(sub_dict['rec_ids']))
            output += "{0}{1}\n".format(" |" * level, ", ".join(sub_dict['rec_ids'])) if print_ids else ""
            output += to_string(sub_dict, level + 1, max_level, print_ids)
    return output


def main():
    import argparse

    parser = argparse.ArgumentParser(prog="split_gb_by_taxa", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("gb_file", nargs="?", default=sys.stdin,
                        help="Path to input file (ensure that taxonomy is annotated)")
    parser.add_argument("-l", "--level", action="store", default=5, type=int,
                        help="Specify how deeply to print output")
    parser.add_argument("-ids", "--print_ids", action="store_true", help="Output the sequence ids")

    in_args = parser.parse_args()

    seqbuddy = Sb.SeqBuddy(in_args.gb_file)
    seqbuddy = Sb.find_repeats(seqbuddy)
    if seqbuddy.repeat_ids:
        print("Error: You have replicate IDs in your file")
        print("------------------------------------------\n")
        for repeat_id in seqbuddy.repeat_ids:
            print(repeat_id)
        sys.exit()

    recs_wo_taxonomy = []
    taxonomy_dict = OrderedDict()
    for rec_id, rec in seqbuddy.to_dict().items():
        if "taxonomy" not in rec.annotations or not rec.annotations['taxonomy']:
            recs_wo_taxonomy.append(rec)
        else:
            if "organism" in rec.annotations:
                rec.annotations['organism'] = re.sub("[ ]+", " ", rec.annotations['organism'])
                if rec.annotations['organism']:
                    organism = rec.annotations['organism'].split()[1:]
                    organism = " ".join(organism)
                    rec.annotations['taxonomy'].append(organism)
            taxonomy_dict = add_tax(rec_id, rec.annotations['taxonomy'], taxonomy_dict)

    print(to_string(taxonomy_dict, level=0, max_level=in_args.level))

if __name__ == '__main__':
    main()
