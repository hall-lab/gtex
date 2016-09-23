#!/usr/bin/env python

import argparse, sys
import re
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
gt_to_cn.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: convert SV genotypes into simple copy number")
    parser.add_argument('-i', '--input', metavar='VCF', dest='vcf_in', type=argparse.FileType('r'), default=None, help='VCF input [stdin]')

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
def gt_to_cn(vcf_in):
    for line in vcf_in:
        if line[:2] == '##':
            print line.rstrip()
            continue

        if line[:6] == "#CHROM":
            # print '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">'
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')

        if not v[2].startswith('LUMPY'):
            print line.rstrip()
            continue

        # parse format tags
        fmt = v[8].split(':')
        ab_idx = -1
        ds_idx = -1
        for i in xrange(len(fmt)):
            if fmt[i] == 'AB':
                ab_idx = i
            elif fmt[i] == 'DS':
                ds_idx = i

        distance_threshold = 1e6
        chrom = v[0]
        pos = int(v[1])
        alt = v[4]
        interchrom = False
        orient = ''
        if "[" in alt or "]" in alt: # BND
            orient = ""
            alt_coord = re.findall('[\[\]]([^\[\]]*)[\[\]]', alt)[0].split(':')
            alt_chrom = alt_coord[0]
            alt_pos = int(alt_coord[1])
            if chrom != alt_chrom:
                interchrom = True
                distance = None
            else:
                interchrom = False
                distance = abs(alt_pos - pos)
                start = pos
                end = alt_pos

            if alt.startswith('['):
                orient = "INV"
            elif alt.startswith(']'):
                if pos > alt_pos:
                    orient = "DEL"
                else:
                    orient = "DUP"
            elif alt.endswith('['):
                if pos > alt_pos:
                    orient = "DUP"
                else:
                    orient = "DEL"
            elif alt.endswith(']'):
                orient = "INV"
            bnd_detail = "BND"
            if interchrom:
                bnd_detail = "INTER_" + bnd_detail
            else:
                if distance > distance_threshold:
                    bnd_detail = "DISTANT_" + bnd_detail + "_" + orient
                else:
                    bnd_detail = "LOCAL_" + bnd_detail + "_" + orient

        if (not interchrom and orient == 'DEL') or v[4].startswith('<DEL'):
            # print 'mod'
            # parse sample tags
            for j in xrange(9,len(v)):
                samp = v[j].split(':')
                samp[ds_idx] = str(-1 * float(samp[ds_idx]))
                v[j] = ':'.join(samp)

        # print the variant line
        print '\t'.join(v)
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    gt_to_cn(args.vcf_in)

    # close the files
    args.vcf_in.close()


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
