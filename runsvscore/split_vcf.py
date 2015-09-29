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
split_vcf.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: clean exclude samples from genome strip file")
    parser.add_argument('-n', '--num_lines', type=int, required=False, default=1000, help='approx. number of lines per file [1000]')
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
def split_vcf(max_lines, vcf_file):

    file_num = 0
    line_counter = 0

    header = []
    # read input VCF
    curr_id = ""
    prev_id = ""

    f = open('split%03d.vcf' % file_num, 'w')

    # read input vcf file
    in_header = True
    while 1:
        line = vcf_file.readline()
        if not line:
            break

        if line[0] == "#":
            header.append(line)
            continue
        elif in_header:
            in_header = False
            for h in header:
                f.write(h)

        # iterate line
        line_counter += 1

        # get the current id
        v = line.rstrip().split('\t')
        prev_id = curr_id
        curr_id = v[2]
        bnd = False
        if 'SVTYPE=BND;' in v[7]:
            bnd = True
            curr_id = curr_id[:-2]

        # write the line to file
        f.write(line)
        
        if line_counter >= max_lines:
            # if the bnd doesn't match the previous record, keep going
            if bnd and curr_id != prev_id:
                continue
            
            # close the previous file
            f.close()

            # open a new file and write the header
            file_num += 1
            f = open('split%03d.vcf' % file_num, 'w')
            for h in header:
                f.write(h)
            line_counter = 0
    f.close()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    split_vcf(args.num_lines, args.input_vcf)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
