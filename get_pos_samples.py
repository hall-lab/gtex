#!/usr/bin/env python

import pysam
import argparse, sys
import math, time, re
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-09-27 09:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
get_pos_samples.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: return a list of positive samples")
    parser.add_argument(metavar='vcf', dest='input_vcf', nargs='?', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = sys.stdin
    # send back the user input
    return args

# primary function
def get_pos_samples(vcf_file):
    vcf_out = sys.stdout

    # read input VCF
    for line in vcf_file:
        if line[0] == '#':
            if line[1] != '#':
                vcf_samples = line.rstrip().split('\t')[9:]
            print line.rstrip()
            continue

        pos_samples = []
        v = line.rstrip().split('\t')
        
        fmt_index_list = v[8].split(':')
        fmt_index = None
        for f in xrange(len(fmt_index_list)):
            if fmt_index_list[f] == 'GT':
                fmt_index = f
                break

        for i in xrange(len(vcf_samples)):
            gt_list = v[i + 9].split(':')
            gt = gt_list[fmt_index]
            if gt != './.' and gt != '0/0':
                pos_samples.append(vcf_samples[i])

        v[7] = v[7] + ';POSSAMP=' + ','.join(pos_samples)

        vcf_out.write('\t'.join(v) + '\n')

    vcf_out.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    get_pos_samples(args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
