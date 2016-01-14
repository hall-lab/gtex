#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import re

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
fuckvcf.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: determine the variant type")
    parser.add_argument('-c', metavar='INT', dest='alt_col', type=int, required=False, default=5, help='ALT column (default 5)')
    parser.add_argument('-d', metavar='INT', dest='distance_threshold', type=int, required=False, default=1000000, help='distance threshold for distant variants [1000000]')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# primary function
def parse_type(alt_col, distance_threshold, infile):
    # distant_threshold = 1000

    for line in infile:
        if line.startswith('#'):
            continue

        svtype = None

        v = line.rstrip().split('\t')
        v_id = v[alt_col - 3]
        v_chrom = v[0]
        v_pos = int(v[1])
        alt = v[alt_col - 1]

        if alt == "<DEL>":
            svtype = "DEL"
        elif alt == "<DUP>":
            svtype = "DUP"
        elif alt.startswith("<DEL:ME:"):
            svtype = "MEI"
        elif alt == "<INV>":
            svtype = "INV"
        elif "[" in alt or "]" in alt: # BND
            orient = ""
            coord = re.findall('[\[\]]([^\[\]]*)[\[\]]', alt)[0].split(':')
            chrom = coord[0]
            pos = int(coord[1])
            if v_chrom != chrom:
                interchrom = True
                distance = None
            else:
                interchrom = False
                distance = abs(pos - v_pos)


            if alt.startswith('['):
                orient = "INV"
            elif alt.startswith(']'):
                if v_pos > pos:
                    orient = "DEL"
                else:
                    orient = "DUP"
            elif alt.endswith('['):
                if v_pos > pos:
                    orient = "DUP"
                else:
                    orient = "DEL"
            elif alt.endswith(']'):
                orient = "INV"
            svtype = "BND"
            if interchrom:
                svtype = "INTER_" + svtype
            else:
                if distance > distance_threshold:
                    svtype = "DISTANT_" + svtype + "_" + orient
                else:
                    svtype = "LOCAL_" + svtype + "_" + orient

        # print line.rstrip()
        # print v_id, svtype, alt

        print  "%s\t%s" % (v_id, svtype)
        # print line.rstrip()
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    parse_type(args.alt_col, args.distance_threshold, args.input)

    # close the input file
    args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
