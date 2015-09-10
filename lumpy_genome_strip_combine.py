#!/usr/bin/env python

import argparse, sys
import gzip
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (colbychiang@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-09-09 15:55 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
lumpy_genome_strip_combine.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    parser.add_argument('-l', '--lumpy_vcf', required=True, type=str, help='LUMPY VCF')
    parser.add_argument('-g', '--gs_vcf', required=True, type=str, help='Genome STRiP VCF')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# primary function
def combine(lumpy_vcf, gs_vcf):
    for line in lumpy_vcf:
        print line.rstrip()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.lumpy_vcf.endswith('.gz'):
        lumpy_vcf = gzip.open(args.lumpy_vcf, 'rb')
    else:
        lumpy_vcf = open(args.lumpy_vcf, 'r')


    if args.gs_vcf.endswith('.gz'):
        gs_vcf = gzip.open(args.gs_vcf, 'rb')
    else:
        gs_vcf = open(args.gs_vcf, 'r')

    # call primary function
    combine(lumpy_vcf, gs_vcf)

    # close the input file
    args.vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
