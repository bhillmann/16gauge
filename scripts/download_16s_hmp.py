#!/usr/bin/env python
import urllib.request
import pandas as pd
import argparse
import sys
import os
import zlib

def make_arg_parser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='If nothing is given, then stdin, else the input file')
    return parser


def stream_gzip_decompress(stream):
    dec = zlib.decompressobj(32 + zlib.MAX_WBITS)  # offset 32 to skip the header
    for chunk in stream:
        rv = dec.decompress(chunk)
        if rv:
            yield rv


def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    with open(args.input) if args.input else sys.stdin as inf:
        df = pd.read_csv(inf)

    file_names = df['SRS_SampleID']
    for file in file_names:
        req = urllib.request.Request('ftp://public-ftp.hmpdacc.org/HM16STR/by_sample/%s.fsa.gz' % file)
        with urllib.request.urlopen(req, 'rb') as ftp_stream:
            with open('%s.fsa' % file, 'wb') as outf:
                outf.writelines(stream_gzip_decompress(ftp_stream))

if __name__ == '__main__':
    main()
