#!/usr/bin/env python
import os
import re
import sys
import argparse

def make_arg_parser():
    parser = argparse.ArgumentParser(description='Aggregate NBCI taxon id counts for a Kraken run.')
    parser.add_argument('-o', '--output', type=str, help='The output file.')
    return parser

def read_fasta(f):
    title = None
    data = None
    for line in f:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield title.strip(), data


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    with open(args.output, 'w') if args.output else sys.stdout as outf:
        for filename in os.listdir(os.getcwd()):
                if filename.endswith('.fsa'):
                    with open(filename) as inf:
                        fasta_itr = read_fasta(inf)
                        for title, data in fasta_itr:
                            m = re.search(r"primer=(.*?) ", title).group(1)
                            if m == 'V3-V5':
                                outf.write('>%s %s\n' % (os.path.splitext(os.path.basename(filename))[0], title))
                                outf.write('%s\n' % data)

if __name__ == '__main__':
    main()