#!/usr/bin/env python

import argparse, sys
import numpy as np
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-07-14 11:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
expr_stats.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: select columns from a file by header names")
    parser.add_argument('-c', '--col', metavar='FILE', required=False, type=argparse.FileType('r'), default=None, help='list of column headers to extract [all]')
    parser.add_argument('-l', '--leading', metavar='INT', required=False, type=int, default=0, help='number of leading columns to print [0]')
    parser.add_argument('-p', '--pass', metavar='STR', dest='pass_prefix', required=False, default=None, help='prefix for comment lines in INPUT to pass unfiltered')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='phenotype file')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

def z_score(value, mean, stdev):
    z = (value - mean) / stdev
    return z

# primary function
def extract_cols(col, lead_cols, pass_prefix, source):
    # get_columns = range(lead_cols)
    if col is not None:
        select = []
        for line in col:
            select.append(line.rstrip())
    
    in_header = True
    for line in source:
        if pass_prefix is not None and  line.startswith(pass_prefix):
            print line.rstrip()
            continue
        v = line.rstrip().split('\t')
        if in_header:
            if col is not None:
                column_map = {c: v.index(c) for c in v}
                get_columns = [column_map[x] for x in select if x in column_map]
                in_header = False
            else:
                get_columns = range(lead_cols, len(v))
                in_header = False
            
            print '\t'.join(v[x] for x in range(lead_cols) + get_columns)
            continue
        data = map(float, [v[x] for x in get_columns])
        # print '\t'.join(v[x] for x in get_columns)
        u = np.mean(data)
        sigma = np.std(data)

        z_list = [z_score(x, u, sigma) for x in data]
        print '\t'.join(map(str, [v[x] for x in xrange(lead_cols)] + z_list))

    source.close()
    if col is not None:
        col.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    extract_cols(args.col, args.leading, args.pass_prefix, args.input)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
