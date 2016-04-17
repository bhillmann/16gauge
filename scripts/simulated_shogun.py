#!/usr/bin/env python
from __future__ import print_function, division

import argparse
from collections import defaultdict
import sys
import pandas as pd
import os
import numpy as np



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

    names = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    with sys.stdin if args.input == '-' else open(args.input, 'r') as inf:
        df = pd.read_csv(inf)
        for i in range(1, 8):

            dd = dict()
            for column in df.columns:
                if column[0] == '#':
                    d = defaultdict(int)
                    for lineage, count in zip(df.ix[:, 0], df[column]):
                        if not np.isnan(count):
                            l = lineage.split(';')
                            if len(l) >= i:
                                d[';'.join(l[:i])] += int(count)
                    dd[column] = d
            df_temp = pd.DataFrame(dd)
            df_temp.to_csv(os.path.join(args.output, names[i-1] + '.csv'))


if __name__ == '__main__':
    main()