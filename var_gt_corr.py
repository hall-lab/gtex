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
    parser.add_argument('-a', '--a_vcf', metavar='VCF', required=True, type=str, default=None, help='VCF')
    parser.add_argument('-b', '--b_vcf', metavar='VCF', required=True, type=str, default=None, help='VCF')
    parser.add_argument('-v', '--variants', dest='variants_file', metavar='FILE', required=True, type=argparse.FileType('r'), help='Variant pairs to compare')
    parser.add_argument('-af', '--a_field', metavar='STR', default='GT', help='specify genotyping format field [GT]')
    parser.add_argument('-bf', '--b_field', metavar='STR', default='GT', help='specify genotyping format field [GT]')


    # parser.add_argument('-b', '--variants', metavar='FILE', dest='variants_file', type=argparse.FileType('r'), default=None, required=False, help='list of variants to include')
    parser.add_argument('-s', '--samples', metavar='FILE', dest='samples_file', type=argparse.FileType('r'), default=None, required=False, help='list of samples to include')

    # # parser.add_argument('-c', '--covar', metavar='FILE', dest='covar', type=argparse.FileType('r'), default=None, required=True, help='tab delimited file of covariates')
    # # parser.add_argument('-v', '--max_var', metavar='FLOAT', dest='max_var', type=float, default=0.1, help='maximum genotype variance explained by covariates for variant to PASS filtering [0.1]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    # if args.vcf_in == None:
    #     if sys.stdin.isatty():
    #         parser.print_help()
    #         exit(1)
    #     else:
    #         args.vcf_in = sys.stdin
    # send back the user input
    return args

# primary function
def var_gt_corr(a_vcf,
                b_vcf,
                a_field,
                b_field,
                a_vars,
                b_vars,
                samp_set):

    X = {} # dict of genotypes for each sample, key is variant id
    # var_ids = []
    samp_cols = []

    for line in a_vcf:
        if line[:2] == '##':
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            for i in xrange(9,len(v)):
                if v[i] in samp_set or len(samp_set) == 0:
                    samp_cols.append(i)
            continue

        var_id = v[2]
        
        if var_id not in a_vars:
            continue

        # read the genotypes
        if a_field == 'GT':
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
                if fmt[i] == a_field:
                    field_idx = i
                    break
            if field_idx == -1:
                sys.stderr.write("Format field '%s' not found for variant %s\n" % (a_field, v[2]))
                exit(1)

            gt_list = []
            for i in samp_cols:
                gt_str = v[i].split(':')[field_idx]

                # if no info for the field, fall back to regular genotype
                if gt_str == '.' or '.' in v[i].split(':')[0]:
                    gt_list.append(-1)
                else:
                    gt_list.append(float(gt_str))

            X[var_id] = gt_list
        print gt_list

    if len(a_vars) != len(X):
        sys.stderr.write("Warning, missing variants\n")
        exit(1)

    # empty array of r values (correlation)
    R = [[0.0] * len(var_list) for i in xrange(len(var_list))]

    for i in xrange(len(var_list)):
        for j in xrange(i,len(var_list)):
            var_pair = np.array([X[var_list[i]], X[var_list[j]]])

            # remove missing genotypes
            var_pair = var_pair[:, var_pair[0]!=-1]

            # ensure non-uniformity in genotype and read depth
            if len(np.unique(var_pair[0,:])) > 1 and len(np.unique(var_pair[1,:])) > 1:
                # calculate regression
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(var_pair)
            else:
                r_value = 'nan'
            R[i][j] = r_value
            R[j][i] = r_value

            # print var_list[i], var_list[j], r_value
    for i in xrange(len(R)):
        print '\t'.join(['%0.6g' % x ** 2 for x in R[i]])

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # get list of variants to examine
    a_vars = []
    b_vars = []
    if args.variants_file is not None:
        for line in args.variants_file:
            v = line.rstrip().split('\t')
            a_vars.append(v[0])
            b_vars.append(v[1])
        args.variants_file.close()

    # get list of samples to examine
    samp_set = set()
    if args.samples_file is not None:
        for line in args.samples_file:
            v = line.rstrip().split('\t')
            samp_set.add(v[0])
        args.samples_file.close()

    # open the VCF files
    if args.a_vcf.endswith('.gz'):
        a_vcf = gzip.open(args.a_vcf, 'rb')
    else:
        a_vcf = open(args.a_vcf, 'r')

    if args.b_vcf.endswith('.gz'):
        b_vcf = gzip.open(args.b_vcf, 'rb')
    else:
        b_vcf = open(args.b_vcf, 'r')

    # call primary function
    var_gt_corr(a_vcf,
                b_vcf,
                args.a_field,
                args.b_field,
                a_vars,
                b_vars,
                samp_set)

    # close the files
    args.a_vcf.close()
    args.b_vcf.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
