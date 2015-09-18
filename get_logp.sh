#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]]
then
    echo "usage: $0 <TISSUE> <EGENE>"
    exit 1
fi

TISSUE=$1
EGENE=$2

zjoin -a tissues/$TISSUE/$EGENE/vars.txt -b <(zcat ../tissues/$TISSUE/${TISSUE}.joint.nom.eqtl.txt.gz | tr ' ' '\t') -1 1 -2 2 | awk -v EGENE=$EGENE '{ if ($2==EGENE) print $1"\t"log($6)/log(10)*-1 }'