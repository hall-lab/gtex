#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
# import numpy as np
# from scipy import stats
import collections
# from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-03-27 09:43 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
template.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: map cadd/svscores to percentiles from a cdf")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='file')
    parser.add_argument('--sv',
                        metavar='FILE', dest='sv_cdf_path',
                        required=True,
                        type=str, default=None,
                        help='SV cdf file')
    parser.add_argument('--snv',
                        metavar='FILE', dest='snv_cdf_path',
                        required=True,
                        type=str, default=None,
                        help='SNV cdf file')
    parser.add_argument('--indel',
                        metavar='FILE', dest='indel_cdf_path',
                        required=True,
                        type=str, default=None,
                        help='indel cdf file')

    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

def ptile(cdf, score):
    for p in reversed(cdf):
        if score >= cdf[p]:
            return p

# primary function
def cdf(input_file, sv_cdf_file, snv_cdf_file, indel_cdf_file):
    sv_ptile = collections.OrderedDict()
    snv_ptile = collections.OrderedDict()
    indel_ptile = collections.OrderedDict()

    for line in sv_cdf_file:
        v = map(float,line.rstrip().split('\t'))
        sv_ptile[v[0]] = v[1]

    for line in snv_cdf_file:
        v = map(float,line.rstrip().split('\t'))
        snv_ptile[v[0]] = v[1]

    for line in indel_cdf_file:
        v = map(float,line.rstrip().split('\t'))
        indel_ptile[v[0]] = v[1]

    for line in input_file:
        v = line.rstrip().split('\t')
        var_type = v[2]
        score = float(v[3])

        if var_type == 'SV':
            cdf = sv_ptile
        elif var_type == 'SNV':
            cdf = sv_ptile
        elif var_type == 'INDEL':
            cdf = indel_ptile

        print '\t'.join(map(str, v + [ ptile(cdf, score) ]))

        


    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # get permutation data
    sv_cdf_file = get_file(args.sv_cdf_path)
    snv_cdf_file = get_file(args.snv_cdf_path)
    indel_cdf_file = get_file(args.indel_cdf_path)

    # call primary function
    cdf(
        input_file,
        sv_cdf_file,
        snv_cdf_file,
        indel_cdf_file
        )

    # close the files
    input_file.close()
    


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
