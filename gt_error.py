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
gt_error.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: inject genotyping error into a fraction of samples")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')
    parser.add_argument('-v', '--variants', metavar='FILE', dest='variants_file', type=argparse.FileType('r'), default=None, required=False, help='list of variants to include')
    parser.add_argument('-f', '--field', metavar='STR', dest='field', default='GT', help='specify genotyping format field [GT]')
    parser.add_argument('-r', '--fraction', metavar='FLOAT', dest='fraction', type=float, default=1.0, help='fraction of genotypes to perturb [1.0]')
    parser.add_argument('-s', '--seed', metavar='INT', type=int, required=False, help='Seed for random number generator')

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
def gt_error(vcf_in, var_list, field, fraction):
    X = {} # dict of genotypes for each sample, key is variant id
    # var_ids = []
    if len(var_list):
        has_var_list = True
    else:
        has_var_list = False

    err_gt_idx = []

    for line in vcf_in:
        if line[:2] == '##':
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')

        if line[0] == "#":
            print line.rstrip()
            num_err_samp = int(round(fraction * (len(v) - 9)))
            if num_err_samp == 0:
                sys.stderr.write('Warning: fraction too low, no genotypes perturbed\n')
            continue

        if has_var_list:
            if v[2] not in var_list:
                continue
        else:
            var_list.append(v[2])

        var_id = v[2]


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

        # get the true genotypes from the input
        gt_list = []
        for i in xrange(9,len(v)):
            gt_str = v[i].split(':')[field_idx]
            gt_list.append(gt_str)

        # inject the error into a new genotype list
        gt_err_list = list(gt_list)
        err_gt_idx = sorted(random.sample(xrange(len(v)-9), num_err_samp))
        for i in err_gt_idx:
            gt_err_list[i] = random.choice(gt_list)

        # insert modified genotypes into the VCF
        for i in err_gt_idx:
            gt_split = v[i + 9].split(':')
            gt_split[field_idx] = gt_err_list[i]
            v[i + 9] = ':'.join(gt_split)

        # print the variant line
        print '\t'.join(v)
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # Seed randomness
    if args.seed != None:
        random.seed(args.seed)
    else:
        random.seed()

    # get list of variants to examine
    var_list = []
    if args.variants_file is not None:
        for line in args.variants_file:
            var_list.append(line.rstrip())
        args.variants_file.close()

    # call primary function
    gt_error(args.vcf_in,
             var_list,
             args.field,
             args.fraction)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
