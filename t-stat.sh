#!/bin/bash


if [[ -z $1 ]]
then
    echo "usage: $0 <LIST>"
    exit 1
fi

# job index
EGENE=`cat ../integrated.eqtls/integrated.eqtls.short.txt | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 2`
TISSUE=`cat ../integrated.eqtls/integrated.eqtls.short.txt | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 3`

# echo -e "$LSB_JOBINDEX\t$EGENE\t$TISSUE"

cat tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt \
    | cut -f 4,5 \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/qnorm.py \
    | paste tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt - \
    | gzip -c \
    > tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt.gz

