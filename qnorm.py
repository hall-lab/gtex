#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
import numpy as np
from scipy import stats
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-09-10 14:53 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
var_gt_corr.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: correlate variants and genotypes")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    return args

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    for line in args.input:

        v = map(float, line.rstrip().split('\t'))
        p_value = v[0]
        beta = v[1]
        
        x = p_value / 2.0

        if beta >= 0:
            z = -stats.norm.ppf(x)
        else:
            z = stats.norm.ppf(x)

        print z

    # close the files
    args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
