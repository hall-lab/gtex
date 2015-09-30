#!/usr/bin/env python

import argparse, sys
import math, time, re
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-09-27 09:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcfToBed.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: convert SV VCF to BED file")
    parser.add_argument(metavar='vcf', dest='input_vcf', nargs='?', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = sys.stdin
    # send back the user input
    return args

# primary function
def vcf_to_bed(vcf_file):
    # read input VCF
    for line in vcf_file:
        if line[0] == '#':
            if line[1] != '#':
                vcf_samples = line.rstrip().split('\t')[9:]
                print '\t'.join(['#CHROM',
                                'POS_START',
                                'POS_END',
                                'ID',
                                'REF',
                                'ALT',
                                'QUAL',
                                'FILTER',
                                'INFO',
                                'FORMAT'] +
                                vcf_samples)
                continue
                                
            print line.rstrip()
            continue

        v = line.rstrip().split('\t')

        info_split = [i.split('=') for i in v[7].split(';')]
        for i in info_split:
            if len(i) == 1:
                i.append(True)
        info = dict(info_split)

        bed_list = []

        if info['SVTYPE'] == 'BND':
            chrom = v[0]
            start = int(v[1]) - 1
            end = int(v[1]) - 1
            event = info['EVENT']
            bed = [chrom, start, end, event]
            bed_list.append(bed)

        elif info['SVTYPE'] == 'INV':
            chrom_1 = v[0]
            start_1 = int(v[1])
            end_1 = int(v[1]) + 1
            event_1 = v[2]
            bed_1 = [chrom_1, start_1, end_1, event_1]

            chrom_2 = v[0]
            start_2 = int(info['END']) - 1
            end_2 = int(info['END'])
            event_2 = v[2]
            bed_2 = [chrom_2, start_2, end_2, event_2]

            bed_list.append(bed_1)
            bed_list.append(bed_2)
            
        else:
            chrom = v[0]
            start = int(v[1])
            end = int(info['END'])
            event = v[2]
            bed = [chrom, start, end, event]
            bed_list.append(bed)

        for b in bed_list:
            print '\t'.join(map(str, b)) + '\t' +  '\t'.join(v[5:])
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    vcf_to_bed(args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
