#!/usr/bin/env python

import argparse, sys
# import math, time, re
import numpy as np
from scipy import stats
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-10-22 11:25 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
scale_gt.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: continuous ld")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-v', '--variants', metavar='FILE', dest='variants_file', type=argparse.FileType('r'), default=None, required=False, help='list of variants to include')
    parser.add_argument('-s', '--samples', metavar='FILE', dest='samples_file', type=argparse.FileType('r'), default=None, required=False, help='list of samples to include')
    parser.add_argument('-f', '--field', metavar='STR', dest='field', default='DS', help='specify genotyping format field [DS]')
    parser.add_argument('-r', '--range', metavar='STR', dest='range', default='0,2', help='range to scale over [0,2]')
    # parser.add_argument('-c', '--covar', metavar='FILE', dest='covar', type=argparse.FileType('r'), default=None, required=True, help='tab delimited file of covariates')
    # parser.add_argument('-v', '--max_var', metavar='FLOAT', dest='max_var', type=float, default=0.1, help='maximum genotype variance explained by covariates for variant to PASS filtering [0.1]')

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
def scale_gt(vcf_in, var_list, samp_set, field, scale_range):
    # X = {} # dict of genotypes for each sample, key is variant id
    # var_ids = []
    samp_cols = []

    for line in vcf_in:
        if line[:2] == '##':
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            for i in xrange(9,len(v)):
                if v[i] in samp_set or len(samp_set) == 0:
                    samp_cols.append(i)
            print line.rstrip()
            continue
        
        if v[2] not in var_list and len(var_list):
            continue

        var_id = v[2]

        # print v[:6]
        # read the genotypes
        fmt = v[8].split(':')
        field_idx = -1
        for i in xrange(len(fmt)):
            if fmt[i] == field:
                field_idx = i
                break
        if field_idx == -1:
            sys.stderr.write("Format field '%s' not found for variant %s\n" % (field, v[2]))
            exit(1)

        fmt_list = []
        for i in samp_cols:
            fmt_str = v[i].split(':')[field_idx]
            fmt_list.append(fmt_str)

        min_fmt = np.min([float(x) for x in fmt_list if x != "."])
        max_fmt = np.max([float(x) for x in fmt_list if x != "."])
        # print 'min_fmt', min_fmt
        # print 'max_fmt', max_fmt
        
        scalar = (max_fmt - min_fmt) / 2.0
        shift = min_fmt
        # print scalar, shift

        # print 'fmt', fmt_list

        fmt_scaled = [(float(x) - shift) / scalar for x in fmt_list]

        # print 'scaled', fmt_scaled

        # print 'min_scaled', np.min([float(x) for x in fmt_scaled if x != "."])
        # print 'max_scaled', np.max([float(x) for x in fmt_scaled if x != "."])

        for i in samp_cols:
            # print 'before', v[i]
            fmt_split = v[i].split(':')
            fmt_split[field_idx] = "%0.6f" % fmt_scaled[i - 9]
            v[i] = ':'.join(fmt_split)
        
        print '\t'.join(v)

        # X[var_id] = fmt_list

    # if len(var_list) != len(X):
    #     sys.stderr.write("Warning, missing variants\n")
    #     exit(1)


    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # get list of variants to examine
    var_list = []
    if args.variants_file is not None:
        for line in args.variants_file:
            var_list.append(line.rstrip())
        args.variants_file.close()

    # get list of samples to examine
    samp_set = set()
    if args.samples_file is not None:
        for line in args.samples_file:
            v = line.rstrip().split('\t')
            samp_set.add(v[0])
        args.samples_file.close()

    scale_range = map(int, args.range.split(','))

    # call primary function
    scale_gt(args.vcf_in, var_list, samp_set, args.field, scale_range)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
