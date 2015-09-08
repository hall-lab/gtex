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
pythonTemplate.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    # parser.add_argument('-x', '--exclude_trunc', required=False, action='store_true', help='excludes truncation scores')
    parser.add_argument('vcf', nargs='?', type=argparse.FileType('r'), default=None, help='VCF file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf = sys.stdin

    # send back the user input
    return args

# primary function
def svscore_consol(vcf):

    svscore_no_trunc_list = ['SVSCOREMAX_SPAN', 'SVSCOREMAX_LEFT', 'SVSCOREMAX_RIGHT']
    svscore_trunc_list = ['SVSCOREMAX_LTRUNC', 'SVSCOREMAX_RTRUNC']

    info_list = []
    for line in vcf:
        if line.startswith("#"):
            if line.startswith("##INFO"):
                a = line.rstrip().split("<")
                b = a[1].split(',')[0][3:]
                info_list.append(b)
            elif line.startswith("#CHROM"):
                if 'SVSCORE' not in info_list:
                    print '##INFO=<ID=SVSCORE,Number=1,Type=Float,Description="Max of all SVScores">'
                if 'SVSCORE_NOTRUNC' not in info_list:
                    print '##INFO=<ID=SVSCORE_NOTRUNC,Number=1,Type=Float,Description="Max of all SVScores excluding truncation">'

            print line.rstrip()

        else:
            svscore_max = 0
            svscore_no_trunc_max = 0
            info = {}
            v = line.rstrip().split('\t')

            for i in v[7].split(';'):
                j = i.split('=')
                if len(j) == 2:
                    info[j[0]] = j[1]
                else:
                    info[j[0]] = True

            # get the max score
            for s in svscore_no_trunc_list:
                try:
                    svscore_max = max(svscore_max, float(info[s]))
                    svscore_no_trunc_max = max(svscore_no_trunc_max, float(info[s]))
                except KeyError:
                    continue

            # get the max score excluding truncation
            for s in svscore_trunc_list:
                try:
                    svscore_max = max(svscore_max, float(info[s]))
                except KeyError:
                    continue

            info['SVSCORE'] = str(svscore_max)
            info['SVSCORE_NOTRUNC'] = str(svscore_no_trunc_max)

            info_string = ';'.join(x + "=" + info[x] for x in info if info[x] != True)
            for x in info:
                if info[x] == True:
                    info_string = info_string + ';' + x
            
            print '\t'.join(v[:7] +
                            [info_string] +
                            v[8:])
            
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    svscore_consol(args.vcf)

    # close the input file
    args.vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
