#!/usr/bin/env python

import argparse, sys
# import math, time, re
import numpy as np
# from scipy import stats
# from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-07-09 11:25 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
grm.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: generate a genetic relatedness matrix from a VCF")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-v', '--variants', metavar='FILE', dest='variants_file', type=argparse.FileType('r'), default=None, required=False, help='list of variants to include')
    parser.add_argument('-s', '--samples', metavar='FILE', dest='samples_file', type=argparse.FileType('r'), default=None, required=False, help='list of samples to include')
    parser.add_argument('-f', '--field', metavar='STR', dest='field', default='GT', help='specify genotyping format field [GT]')
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
def make_grm(vcf_in, var_set, samp_set, field):
    X = [] # matrix of genotypes for each sample
    var_ids = []
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
        
        if v[2] not in var_set:
            continue

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

            X.append(gt_list)
            var_ids.append(v[2])
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
                if gt_str == '.' or '.' in v[i].split(':')[0]:
                    gt_list.append(-1)
                else:
                    gt_list.append(float(gt_str))

            X.append(gt_list)
            var_ids.append(v[2])

    N = len(X) # number of variants
    S = len(X[0]) # number of samples

    p = [] # population allele frequency of alternate allele
    d = [] # denominator to normalize variant
    for i in xrange(len(X)): # each i is a different variant
        gt = X[i]
        informative_gt = [g for g in X[i] if g != -1]
        p_i = sum(informative_gt) / (2.0 * len(informative_gt))
        p.append(p_i)

        # fill missing genotypes with the population mean
        if X[i] == -1: X[i] = p_i

        diff = [(g - p_i) for g in gt]
        d_i = sum([d_j ** 2 for d_j in diff]) ** 0.5
        d.append(d_i)
                         
    for j in xrange(S):
        for k in xrange(j + 1):
            gr = 0.0
            num_obs = 0
            # print j,k
            # print grm[j][k]
            for i in xrange(N):
                if (p[i] > 0 and p[i] < 1
                    and  X[i][j] != -1 and X[i][k] != -1):
                    num_obs += 1
                    gr += (X[i][j] - 2 * p[i]) * (X[i][k] - 2 * p[i]) / ( 2 * p[i] * (1 - p[i]))
            gr = gr / float(N)

            print "%s\t%s\t%s\t%.6g" % (j + 1, k + 1, num_obs, gr)

    # for i in xrange(len(grm)):
    #     print '\t'.join(map(str, grm[i]))
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # get list of variants to examine
    var_set = set()
    if args.variants_file is not None:
        for line in args.variants_file:
            var_set.add(line.rstrip())
        args.variants_file.close()

    # get list of samples to examine
    samp_set = set()
    if args.samples_file is not None:
        for line in args.samples_file:
            v = line.rstrip().split('\t')
            samp_set.add(v[0])
        args.samples_file.close()

    # call primary function
    make_grm(args.vcf_in, var_set, samp_set, args.field)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
