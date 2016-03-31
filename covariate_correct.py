#!/usr/bin/env python

import argparse, sys
import gzip
from argparse import RawTextHelpFormatter
from scipy import stats
import numpy as np
import statsmodels.api as sm
# from sklearn import datasets, linear_model

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-03-30 18:18 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
covariate_correct.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: correct a matrix with residualization by linear regression of covariates")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None, required=False,
                        help='matrix of uncorrected values')
    parser.add_argument('-c', '--covar',
                        metavar='FILE', dest='covar_path',
                        required=True,
                        type=str, default=None,
                        help='covariates')
    parser.add_argument('-C', '--skip_cols',
                        metavar='INT', dest='skip_cols',
                        required=False,
                        type=int, default=0,
                        help='number of leading columns to skip [0]')
    parser.add_argument('-z', '--z_norm',
                        dest='z_norm',
                        action='store_true',
                        help = 'report row-normalized z-scores')
    # parser.add_argument('-R', '--skip_rows',
    #                     metavar='INT', dest='skip_rows',
    #                     required=False,
    #                     type=int, default=0,
    #                     help='number of leading rows to skip [0]')


    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

# primary function
def resid_lin_reg(
    matrix_file,
    covar_file,
    skip_cols,
    z_norm
    ):

    # read through the covariates
    in_header = True
    cov_header = []
    cov_list = []
    cov_list_sorted = []
    for line in covar_file:
        v = line.rstrip().split('\t')
        if in_header:
            cov_header = v[1:]
            in_header = False
            continue
        cov_list.append(map(float, v[1:]))

    # read through the matrix
    in_header = True
    mat_header = []
    for line in matrix_file:
        v = line.rstrip().split('\t')
        if in_header:
            print line.rstrip()
            mat_header = v[skip_cols:]
            in_header = False
            # rearrange the covariate table so it matches the matrix
            reorder = [cov_header.index(x) for x in mat_header]
            for row in cov_list:
                cov_list_sorted.append([row[x] for x in reorder])

            # X is numpy array of sorted covariates
            X = np.transpose(np.array(cov_list_sorted))
            # print X
            continue

        # expression
        y = np.array(map(float, v[skip_cols:]))
        results = sm.OLS(y, X).fit()

        # Inspect the results
        # print results.summary()

        # store residuals of matrix values
        resid = results.resid.tolist()

        # normalize by z-score if requested
        if z_norm:
            corrected = stats.zscore(resid)
        else:
            corrected = resid

        # print results
        print '\t'.join(v[:skip_cols] + map(str, corrected))

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # get permutation data
    covar_file = get_file(args.covar_path)

    # call primary function
    resid_lin_reg(
        input_file,
        covar_file,
        args.skip_cols,
        args.z_norm
        )

    # close the files
    input_file.close()
    covar_file.close()
    


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
