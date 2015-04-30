#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
pheno_samples.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: extract samples in VCF from phenotype file")
    parser.add_argument('-v', '--vcf', required=True, type=argparse.FileType('r'), help='VCF containing desired samples')
    parser.add_argument('-c', '--col', required=False, type=int, default=0, help='number of leading columns to print [0]')
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

# primary function
def pheno_samples(vcf, lead_cols, pheno):
    get_columns = range(lead_cols)
    for line in vcf:
        if line[:2]=="##":
            continue
        elif line[0] == "#":
            get_samples = line.rstrip().split('\t')[9:]
            break
    
    in_header = True
    for line in pheno:
        v = line.rstrip().split('\t')
        if in_header:
            for i in xrange(lead_cols, len(v)):
                if v[i] in get_samples:
                    get_columns.append(i)
            in_header = False
        print '\t'.join(v[x] for x in get_columns)

    # in_header = False
    # print get_columns

    pheno.close()
    vcf.close()    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    pheno_samples(args.vcf, args.col, args.input)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
