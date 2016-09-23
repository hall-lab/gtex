#!/usr/bin/env python

import argparse, sys
import random
# import math, time, re
# import numpy as np
# from scipy import stats
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-09-09 10:25 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
gt_to_cn.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: convert SV genotypes into simple copy number")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf_in == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf_in = sys.stdin
    # send back the user input
    return args

# primary function
def gt_to_cn(vcf_in):
    for line in vcf_in:
        if line[:2] == '##':
            print line.rstrip()
            continue

        if line[:6] == "#CHROM":
            print '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">'
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')

        v[8] = v[8] + ":DS"

        # parse the alt
        is_del = False
        is_cn = False
        alt = v[4].split(',')
        for j in xrange(len(alt)):
            if alt[j].startswith('<CN'):
                is_cn = True
                alt[j] = int(alt[j].strip('<CN>'))

        if not is_cn:
            if 'SVTYPE=DEL' in v[7]:
                is_del = True
            alt = [1, 0]

        # tack on the hom ref genotype
        if is_cn:
            alt = [1] + alt



        # parse the sample GTs
        for i in xrange(9,len(v)):
            out_gt = 0
            gt_str = v[i].split(':')[0]
            if '.' in gt_str:
                out_gt = '.'
                continue
            else:
                sep = '/'
                if sep not in gt_str:
                    sep = '|'

                if is_cn or is_del:
                    out_gt = 0
                    for g in map(int, gt_str.split(sep)):
                        out_gt += alt[g]
                else:
                    out_gt = sum(map(int, gt_str.split(sep)))

            # append to the sample column
            v[i] = v[i] + ':' + str(out_gt)

        # print the variant line
        print '\t'.join(v)
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    gt_to_cn(args.vcf_in)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
