#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
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
    parser.add_argument('-a', '--algorithm', metavar='STR', dest='algorithm', default='mott', help='algorithm to use (mott, visscher) [mott]')
    parser.add_argument('-z', '--znorm', dest='znorm', action='store_true', help='z-normalize genotypes prior to GRM')
    parser.add_argument('-o', '--out', metavar='STR', dest='out_prefix', required=True, help='output file prefix')
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

# mott algorithm
def mott(X, N, p, j, k):
    gr = 0.0
    a = 0.0
    b = 0.0
    c = 0.0
    num_obs = 0

    for i in xrange(N):
        if (len(set(X[i])) != 1
            and  X[i][j] != -1 and X[i][k] != -1):
            num_obs += 1

            a += (X[i][j] - 2 * p[i]) * (X[i][k] - 2 * p[i])
            b += (X[i][j] - 2 * p[i]) ** 2
            c += (X[i][k] - 2 * p[i]) ** 2
    gr = a / ((b * c) ** 0.5)
    return (gr, num_obs)

# visscher algorithm
def visscher(X, N, p, j, k):
    gr = 0.0
    num_obs = 0

    for i in xrange(N):
        if (len(set(X[i])) != 1
            and  X[i][j] != -1 and X[i][k] != -1
            and p[i] > 0 and p[i] < 1):
            num_obs += 1

            gr += (X[i][j] - 2 * p[i]) * (X[i][k] - 2 * p[i]) / ( 2 * p[i] * (1 - p[i]))

    gr = gr / float(N)
    return (gr, num_obs)


# primary function
def make_grm(vcf_in,
             var_set,
             samp_set,
             field,
             algorithm,
             znorm,
             out_prefix):
    out_grm = gzip.open("%s.grm.gz" % out_prefix, 'wb')
    out_id = open("%s.grm.id" % out_prefix, 'w')


    X = [] # matrix of genotypes for each sample
    var_ids = []
    samp_cols = []

    sys.stderr.write("Reading genotypes... ")
    for line in vcf_in:
        if line[:2] == '##':
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            for i in xrange(9,len(v)):
                if v[i] in samp_set or len(samp_set) == 0:
                    samp_cols.append(i)
                    out_id.write("%s\t%s\n" % (v[i], v[i]))
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
                if gt_str == '.':
                    gt_list.append(-1)
                else:
                    gt_list.append(float(gt_str))

            if znorm:
                gt_mean = np.mean(gt_list)
                gt_std = np.std(gt_list)
                if gt_std == 0:
                    gt_list = [0 for gt in gt_list]
                else:
                    gt_list = [(gt - gt_mean) / gt_std for gt in gt_list]

            X.append(gt_list)
            var_ids.append(v[2])
    # close the id file
    out_id.close()

    # done reading genotypes
    sys.stderr.write("done\n")

    sys.stderr.write("Calculating variant statistics... ")
    N = len(X) # number of variants
    S = len(X[0]) # number of samples
    p = [] # population allele frequency of alternate allele
    d = [] # denominator to normalize variant
    for i in xrange(len(X)): # each i is a different variant
        gt = X[i]
        informative_gt = [g for g in X[i] if g != -1]
        try:
            p_i = sum(informative_gt) / (2.0 * len(informative_gt))
        except ZeroDivisionError:
            p_i = -1
        p.append(p_i)

        # fill missing genotypes with the population mean
        if X[i] == -1: X[i] = p_i

        diff = [(g - p_i) for g in gt]
        d_i = sum([d_j ** 2 for d_j in diff]) ** 0.5
        d.append(d_i)
    sys.stderr.write("done\n")

    sys.stderr.write("Calculating genetic relatedness...\n")
    for j in xrange(S):
        sys.stderr.write("%s\n" % (j + 1))
        for k in xrange(j + 1):
            # print j,k
            # print grm[j][k]

            if algorithm == 'mott':
                (gr, num_obs) = mott(X, N, p, j, k)
            elif algorithm == 'visscher':
                (gr, num_obs) = visscher(X, N, p, j, k)

            # print "%s\t%s\t%s\t%.6g" % (j + 1, k + 1, num_obs, gr)
            out_grm.write("%s\t%s\t%s\t%.8g\n" % (j + 1, k + 1, num_obs, gr))

    # close the grm file
    out_grm.close()
    
    # done calculating genetic relatedness
    sys.stderr.write("done\n")

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

    allowed_algorithms = ['mott', 'visscher']
    if args.algorithm not in allowed_algorithms:
        sys.stderr.write('Error: algorithm "%s" not recognized. Choose from [%s]\n' % (args.algorithm, ','.join(allowed_algorithms)))
        exit(1)

    # call primary function
    make_grm(args.vcf_in,
             var_set, samp_set,
             args.field,
             args.algorithm,
             args.znorm,
             args.out_prefix)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
