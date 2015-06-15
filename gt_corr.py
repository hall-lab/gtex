#!/usr/bin/env python

import argparse, sys, copy, gzip
import math, time, re
import numpy
from scipy import stats
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.2 $"
__date__ = "$Date: 2014-04-28 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
gt_corr.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: derive correlation between variant genotypes")
    parser.add_argument('-a', '--vcf_a', metavar='VCF', dest='vcf_a', type=str, required=True, help='VCF input')
    parser.add_argument('-b', '--vcf_b', metavar='VCF', dest='vcf_b', type=str, required=True, help='VCF input b')

    # parse the arguments
    args = parser.parse_args()

    return args


# test whether variant has read depth support
def has_depth_support(var):
    slope_threshold = 0.1
    rsquared_threshold = 0.1
    
    if 'CN' in var.active_formats:
        gt_list = []
        for s in var.sample_list:
            gt_str = var.genotype(s).get_format('GT')
            if '.' in gt_str:
                gt_list.append(-1)
                continue

            sep = '/'
            if sep not in gt_str:
                sep = '|'
            gt_list.append(sum(map(int, gt_str.split(sep))))

        rd_list = map(float, [var.genotype(s).get_format('CN') for s in var.sample_list])
        rd = numpy.array([gt_list, rd_list])

        # remove missing genotypes
        rd = rd[:, rd[0]!=-1]

        # ensure non-uniformity in genotype and read depth
        if len(numpy.unique(rd[0,:])) > 1 and len(numpy.unique(rd[1,:])) > 1:
            # calculate regression
            (slope, intercept, r_value, p_value, std_err) = stats.linregress(rd)
            # print slope, intercept, r_value, var.info['SVTYPE'], var.var_id

            # # write the scatterplot to a file
            # f = open('data/%s_%s_%sbp.txt' % (var.info['SVTYPE'], var.var_id, var.info['SVLEN']), 'w')
            # numpy.savetxt(f, numpy.transpose(rd), delimiter='\t')
            # f.close()
            
            if r_value ** 2 < rsquared_threshold:
                return False

            if var.info['SVTYPE'] == 'DEL':
                slope = -slope

            if slope < slope_threshold:
                return False

            return True
    return False

# primary function
def gt_correlate(vcf_a, vcf_b):
    max_dist = 1e5

    a_vcf_list = []
    b_vcf_list = []

    for a_line in vcf_a:
        if a_line.startswith('#'):
            continue

        a_v = a_line.rstrip().split('\t')
        a_vcf_list.append(a_v)

    for b_line in vcf_b:
        if b_line.startswith('#'):
            continue
        b_v = b_line.rstrip().split('\t')
        b_vcf_list.append(b_v)

    for a_v in a_vcf_list:
        a_id = a_v[2]
        a_chrom = a_v[0]
        a_pos = int(a_v[1])

        a_gt = []
        for i in xrange(9,len(a_v)):
            gt_str = a_v[i].split(':')[0]
            if '.' in gt_str:
                a_gt.append(-1)
                continue
            sep = '/'
            if sep not in gt_str:
                sep = '|'
            a_gt.append(sum(map(int, gt_str.split(sep))))

        # read vcf_b
        for b_v in b_vcf_list:
            b_id = b_v[2]
            b_chrom = b_v[0]
            b_pos = int(b_v[1])

            if (b_chrom != a_chrom
                or abs(b_pos - a_pos > max_dist)):
                continue

            b_gt = []
            for i in xrange(9,len(b_v)):
                gt_str = b_v[i].split(':')[0]
                if '.' in gt_str:
                    b_gt.append(-1)
                    continue
                sep = '/'
                if sep not in gt_str:
                    sep = '|'
                b_gt.append(sum(map(int, gt_str.split(sep))))
            

            glink = numpy.array([a_gt, b_gt])

            # remove missing genotypes
            glink = glink[:, glink[0]!=-1]

            # ensure non-uniformity in genotype and read depth
            if len(numpy.unique(glink[0,:])) > 1 and len(numpy.unique(glink[1,:])) > 1:
                # calculate regression
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(glink)
                print '\t'.join(map(str, [a_id, b_id, r_value ** 2]))

                # # write the scatterplot to a file
                # f = open('data/%s_%s_%sbp.txt' % (var.info['SVTYPE'], var.var_id, var.info['SVLEN']), 'w')
                # numpy.savetxt(f, numpy.transpose(glink), delimiter='\t')
                # f.close()


                # if r_value ** 2 < rsquared_threshold:
                #     return False

                # if var.info['SVTYPE'] == 'DEL':
                #     slope = -slope

                # if slope < slope_threshold:
                #     return False

                # return True

    return


# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.vcf_a.endswith('.gz'):
        vcf_a = gzip.open(args.vcf_a, 'rb')
    else:
        vcf_a = open(args.vcf_a, 'r')

    if args.vcf_b.endswith('.gz'):
        vcf_b = gzip.open(args.vcf_b, 'rb')
    else:
        vcf_b = open(args.vcf_b, 'r')

    # call primary function
    gt_correlate(vcf_a, vcf_b)

    # close the files
    vcf_a.close()
    vcf_b.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
