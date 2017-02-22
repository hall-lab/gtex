#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import math
import numpy as np

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "0.0.2"
__date__ = "$Date: 2015-04-21 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
quantile_norm.py " + __version__ + "\n\
author: " + __author__ + "\n\
description: upper quartile normalize expression data")
    # parser.add_argument('-a', '--abs', action='store_true', help='take absolute values of input before calculating stats')
    # parser.add_argument('-q', '--quantile', type=str, default=0.75, help='quantile for normalization [0.75]')
    parser.add_argument('-r', '--skip_rows', type=int, default=0, help='number of rows to skip [0]')
    parser.add_argument('-c', '--skip_columns', type=int, default=0, help='number of columns to skip [0]')
    parser.add_argument('data', nargs='?', type=argparse.FileType('r'), default=None, help='input data [stdin]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.data == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.data = sys.stdin
    return args

# primary function
def qnorm(data,
          skip_rows,
          skip_columns):

    # store the row and column headers so we can
    # output them at the end
    d_str = np.genfromtxt(data, dtype=str, comments=None)
    row_head = d_str[:skip_rows,:]
    col_head = d_str[:,:skip_columns]

    d = d_str[skip_rows:,skip_columns:].astype(float)
    num_rows = d.shape[0]
    num_cols = d.shape[1]

    d_norm = np.zeros_like(d)
    sorted_sums = np.zeros(num_rows)

    # get the means of the sorted rows
    for j in xrange(num_cols):
        sorted_sums += np.sort(d[:,j])
    sorted_means = (sorted_sums / num_cols).tolist()

    # get the ranked mean value of each sample's expression
    for j in xrange(num_cols):
        index = np.searchsorted(np.sort(d[:,j]), d[:,j])
        d_norm[:,j] = [sorted_means[i] for i in index]
    
    # print the header
    for i in xrange(row_head.shape[0]):
        print '\t'.join(row_head[i,:])

    # print the name columns plus the actual data
    for i in xrange(num_rows):
        print '\t'.join(col_head[i + skip_rows,:].tolist() +
                        ["%0.10f" % x for x in d_norm[i,:]]
                        )
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    qnorm(args.data,
          args.skip_rows,
          args.skip_columns)
    
    # close the file
    args.data.close()

# initialize the script
if __name__ == '__main__':
    try:
        main()
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
