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
    parser.add_argument('-I', '--index', metavar='STR', dest='index', required=False, type=str, help="get LD of each variant against a single index variant")
    parser.add_argument('-l', '--labels', dest='labels', required=False, action='store_true', help='attach labels to LD matrix')
    parser.add_argument('-c', '-columns', dest='columns', required=False, action='store_true', help='display output in column (rather than matrix) format')
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
def ld_continuous(vcf_in, var_list, samp_set, field, alg, index_var, labels, columns):
    X = {} # dict of genotypes for each sample, key is variant id
    # var_ids = []
    samp_cols = []
    if len(var_list):
        has_var_list = True
    else:
        has_var_list = False

    for line in vcf_in:
        if line[:2] == '##':
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            for i in xrange(9,len(v)):
                if v[i] in samp_set or len(samp_set) == 0:
                    samp_cols.append(i)
            continue
        
        if v[2] not in var_list and has_var_list:
            continue
        else:
            var_list.append(v[2])        

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

    if len(var_list)==0:
        var_list = X.keys()

    if len(var_list) != len(X):
        sys.stderr.write("Warning, missing variants\n")
        exit(1)

    if index_var is None:
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

        # print output
        # in column format
        if columns:
            for i in xrange(len(R)):
                for j in xrange(len(R)):
                    if alg == 'r':
                        ld = R[i][j]
                    elif alg == 'r2':
                        ld = R[i][j] **2
                    print '\t'.join(map(str, (var_list[i], var_list[j], ld)))

        # in matrix format
        else:
            if labels:
                print '\t' + '\t'.join(var_list)
            if alg == 'r':
                for i in xrange(len(R)):
                    if labels:
                        sys.stdout.write(var_list[i] + '\t')
                    print '\t'.join(['%0.6g' % x for x in R[i]])
            elif alg == 'r2':
                for i in xrange(len(R)):
                    if labels:
                        sys.stdout.write(var_list[i] + '\t')
                    print '\t'.join(['%0.6g' % x ** 2 for x in R[i]])

    # test against a single variant
    else:
        R_index_var = [None] * len(var_list)
        # for i in var_list.index(index_var):
        i =  var_list.index(index_var)
        for j in xrange(len(var_list)):
            # extract the variant pair from the dictionary
            var_pair = np.array([X[var_list[i]], X[var_list[j]]])

            # calculate regression
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(var_pair)

            R_index_var[j] = r_value

        # print output
        for j in xrange(len(R_index_var)):
            if alg == 'r':
                value = R_index_var[j]
            elif alg == 'r2':
                value = R_index_var[j] ** 2

            print "%s\t%s\t%0.6g" % (var_list[j], index_var, value)

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
    ld_continuous(args.vcf_in, var_list, samp_set, args.field, args.alg, args.index, args.labels, args.columns)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
