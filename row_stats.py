#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import numpy as np

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-08-18 11:38 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
row_stats.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: select columns from a file by header names")
    parser.add_argument('-l', '--leading', metavar='INT', required=False, type=int, default=0, help='number of leading columns to print [0]')
    parser.add_argument('-p', '--pass', metavar='STR', dest='pass_prefix', required=False, default=None, help='prefix for comment lines in INPUT to pass unfiltered')
    parser.add_argument('-s', '--stats', metavar='STR', dest='query_stats', required=True, help='list of stats (mean,median,mode)')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='tab-delimited file')

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

# primary function
def row_stats(lead_cols, pass_prefix, query_stats, source):
    for line in source:
        data = []
        stats = []

        if pass_prefix is not None and  line.startswith(pass_prefix):
            print line.rstrip()
            continue
        else:
            v = line.rstrip().split('\t')
            for i in xrange(lead_cols, len(v)):
                # print v[i]
                try:
                    data.append(float(v[i]))
                except ValueError:
                    continue

            if len(data) == 0:
                for q in query_stats:
                    if q == 'count':
                        s = len(data)
                    else:
                        s = 'NA'
            else:
                for q in query_stats:
                    if q == 'mean':
                        s = np.mean(data)
                    elif q == 'median':
                        s = np.median(data)
                    elif q == 'mode':
                        s = np.mode(data)
                    elif q == 'min':
                        s = min(data)
                    elif q == 'max':
                        s = max(data)
                    elif q == 'sum':
                        s = np.sum(data)
                    elif q == 'product':
                        s = np.prod(data)
                    elif q == 'count':
                        s = len(data)
                    elif q == 'median_col':
                        median = np.median(data)
                        s = data.index(median)
                    stats.append(s)
        
            print '\t'.join(v[x] for x in xrange(lead_cols)) + '\t' + '\t'.join(map(str, stats))

    source.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    query_stats = args.query_stats.split(',')

    allowed_stats = ['mean', 'median', 'min', 'max', 'sum', 'product', 'count', 'median_col']
    for q in query_stats:
        if q not in allowed_stats:
            sys.stderr.write('Error: %s not in allowed stats (%s)\n' % (q, ','.join(allowed_stats)))
            exit(1)

    # call primary function
    row_stats(args.leading, args.pass_prefix, query_stats, args.input)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
