#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
make_pheno_bed.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    # parser.add_argument('-i', '--', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('-m', '--map', 
                        metavar="FILE",
                        dest='sample_map_file',
                        type=argparse.FileType('r'),
                        required=True,
                        help='tab delimited file mapping VCF samples (column 1) to phenotype samples (column 2).\nexample: NA12878\tNA12878_liver')
    parser.add_argument('-i', '--input',
                        metavar="FILE",
                        required=False,
                        type=argparse.FileType('r'),
                        default=None,
                        help='file to read. [stdin]')

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

# generates the phenotype file
def make_pheno_bed(input, sample_map_file):
    sample_map = dict() # map of expression sample name to vcf sample name
    phen_col = dict() # map from vcf sample name to phenotype column number
    in_head = True

    # process sample map
    for line in sample_map_file:
        v = line.rstrip().split("\t")
        sample_map[v[1]] = v[0]
    sample_map_file.close()

    # process the phenotype file
    for line in input:
        v = line.rstrip().split('\t')

        if in_head:
            for i in xrange(4,len(v)):
                if v[i] in sample_map:
                    phen_col[v[i]] = i
            print '\t'.join(['#CHROM', 'START', 'END', 'ID'] +
                            [sample_map[s] for s in sample_map])
            in_head = False
            continue

        chrom = v[2]
        start = str(int(v[3])-1)
        end = str(int(v[3]))
        transcript_id = v[0]

        # output BED line
        print '\t'.join([chrom, start, end, transcript_id] +
                        [v[phen_col[s]] for s in sample_map])

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    make_pheno_bed(args.input, args.sample_map_file)

    # close the input file
    args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
