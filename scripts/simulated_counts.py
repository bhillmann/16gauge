#!/usr/bin/env python
from __future__ import print_function, division

import argparse
from collections import defaultdict
import sys
import pandas as pd
import os

from shogun.taxonomy.ncbi_tree import NCBITree


def make_arg_parser():
    parser = argparse.ArgumentParser(description='Get least common ancestor for alignments in unsorted BAM/SAM file')
    parser.add_argument('-i', '--input', help='The folder containing the SAM files to process.', required=True, type=str)
    parser.add_argument('-o', '--output', help='If nothing is given, then STDOUT, else write to file')
    return parser


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    tree = NCBITree()

    names = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    names_check = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    with sys.stdin if args.input == '-' else open(args.input, 'r') as inf:
        df = pd.read_csv(inf)

        def check(lin):
            split = lin.split(';')
            for dd in zip(split, names_check):
                if not dd[0][0] == dd[1][0]:
                    return False
            return True

        mp_lineages = [tree.mp_lineage(ncbi_taxon_id) for ncbi_taxon_id in df.ncbi_taxon_id]

        mp_lineages_index = [check(lin) for lin in mp_lineages]
        for i in range(1, 8):
            d = defaultdict(int)
            for lineage, count, ind in zip(mp_lineages, df['count'], range(len(mp_lineages))):
                if mp_lineages_index[ind]:
                    d[';'.join(lineage.split(';')[:i])] += count
            dd = {'#40000000': d}
            series = pd.DataFrame(dd)
            series.to_csv(os.path.join(args.output, names[i-1] + '.csv'))


if __name__ == '__main__':
    main()