#!/usr/bin/env python

import argparse, sys
# import math, time, re
import numpy as np
from scipy import stats
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-07-09 11:25 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
ld_continuous.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: continuous ld")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-v', '--variants', metavar='FILE', dest='variants_file', type=argparse.FileType('r'), default=None, required=False, help='list of variants to include')
    parser.add_argument('-s', '--samples', metavar='FILE', dest='samples_file', type=argparse.FileType('r'), default=None, required=False, help='list of samples to include')
    parser.add_argument('-f', '--field', metavar='STR', dest='field', default='GT', help='specify genotyping format field [GT]')
    parser.add_argument('-a', '--alg', metavar='STR', dest='alg', required=True, type=str, help="LD algorithm ('r', 'r2')")
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
def ld_continuous(vcf_in, var_list, samp_set, field, alg):
    X = {} # dict of genotypes for each sample, key is variant id
    # var_ids = []
    samp_cols = []

    for line in vcf_in:
        if line[:2] == '##':
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            for i in xrange(9,len(v)):
                if v[i] in samp_set or len(samp_set) == 0:
                    samp_cols.append(i)
            continue
        
        if v[2] not in var_list and len(var_list):
            continue

        var_id = v[2]

        # read the genotypes
        if field == 'GT':
            gt_list = []
            for i in samp_cols:
                gt_str = v[i].split(':')[0]
                if '.' in gt_str:
                    gt_list.append(-1)
                    continue

                sep = '/'
                if sep not in gt_str:
                    sep = '|'
                gt_list.append(sum(map(int, gt_str.split(sep))))

            X[var_id] = gt_list
        else:
            fmt = v[8].split(':')
            field_idx = -1
            for i in xrange(len(fmt)):
                if fmt[i] == field:
                    field_idx = i
                    break
            if field_idx == -1:
                sys.stderr.write("Format field '%s' not found for variant %s\n" % (field, v[2]))
                exit(1)

            gt_list = []
            for i in samp_cols:
                gt_str = v[i].split(':')[field_idx]

                # if no info for the field, fall back to regular genotype
                if gt_str == '.':
                    gt_list.append(-1)
                else:
                    gt_list.append(float(gt_str))

            X[var_id] = gt_list

    if len(var_list) != len(X):
        sys.stderr.write("Warning, missing variants\n")
        exit(1)

    # empty array of r values (correlation)
    R = [[0.0] * len(var_list) for i in xrange(len(var_list))]

    for i in xrange(len(var_list)):
        for j in xrange(i,len(var_list)):
            # extract the variant pair from the dictionary
            var_pair = np.array([X[var_list[i]], X[var_list[j]]])

            # print var_pair

            # calculate regression
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(var_pair)

            # print 'r_value:', r_value

            R[i][j] = r_value
            R[j][i] = r_value

            # print var_list[i], var_list[j], r_value
    if alg == 'r':
        for i in xrange(len(R)):
            print '\t'.join(['%0.6g' % x for x in R[i]])
    elif alg == 'r2':
        for i in xrange(len(R)):
            print '\t'.join(['%0.6g' % x ** 2 for x in R[i]])


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

    # parse algorithm
    if args.alg not in ('r', 'r2'):
        sys.stderr.write("\nError: algorithm '%s' not supported. Must be 'r' or 'r2'\n\n" % args.alg)
        exit(1)

    # call primary function
    ld_continuous(args.vcf_in, var_list, samp_set, args.field, args.alg)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
