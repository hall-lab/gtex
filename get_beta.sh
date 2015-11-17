#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]] || [[ -z $3 ]]
then
    echo "usage: $0 <TISSUE> <EGENE> <EQTL FILE>"
    exit 1
fi

TISSUE=$1
EGENE=$2
FILE=$3

zjoin -a tissues/$TISSUE/$EGENE/vars.txt -b <(zcat $FILE | tr ' ' '\t') -1 1 -2 2 | awk -v EGENE=$EGENE '{ if ($2==EGENE) print $1"\t"$5 }'