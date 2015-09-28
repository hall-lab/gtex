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
clean_gs.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: clean exclude samples from genome strip file")
    parser.add_argument('-e', '--exclude', type=argparse.FileType('r'), required=True, help='samples to exclude private variants')
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
def clean_gs(exclude_file, vcf_file):
    vcf_out = sys.stdout

    exclude = []
    for line in exclude_file:
        exclude.append(line.rstrip())

    # read input VCF
    for line in vcf_file:
        if line[0] == '#':
            if line[1] != '#':
                vcf_samples = line.rstrip().split('\t')[9:]
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')
        
        fmt_index_list = v[8].split(':')
        fmt_index = None
        for f in xrange(len(fmt_index_list)):
            if fmt_index_list[f] == 'CN':
                fmt_index = f
                break

        # samples
        exclude_gts = []
        include_gts = []

        for i in xrange(len(vcf_samples)):
            gt_list = v[i + 9].split(':')
            cn = gt_list[fmt_index]
            if cn == '.':
                continue
            else:
                cn = int(cn)

            if vcf_samples[i] in exclude:
                exclude_gts.append(cn)
            else:
                include_gts.append(cn)

        if len(set(include_gts)) <= 1:
            v[6] = "RD"

        print '\t'.join(v)

        # print 'excl', exclude_gts
        # print include_gts
        # print vcf.sample_list
        # for sample in vcf.sample_list:
        #     if sample in exclude:
        #         exclude_gts.append(var.genotype(sample).get_format('CN'))
        #     else:
        #         include_gts.append(var.genotype(sample).get_format('CN'))

        # print exclude_gts
        # print include_gts

    #     vcf_out.write(var.get_var_string() + '\n')
    # vcf_out.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    clean_gs(args.exclude, args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
