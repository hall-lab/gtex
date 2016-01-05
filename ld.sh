#!/bin/bash


if [[ -z $1 ]]
then
    echo "usage: $0 <LIST>"
    exit 1
fi

# job index
EGENE=`cat ../integrated.eqtls/integrated.eqtls.short.txt | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 2`
TISSUE=`cat ../integrated.eqtls/integrated.eqtls.short.txt | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 3`

REGION=`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz | awk -v EGENE=$EGENE -v SLOP=10000000 '{ if (\$1==EGENE) { POS_START=\$3-SLOP; POS_END=\$4+SLOP; if (POS_START<0) POS_START=0; print \$2":"POS_START"-"POS_END } }'`;

# echo -e "$LSB_JOBINDEX\t$EGENE\t$TISSUE"

tabix -h ../../gtex.joint_sv_gatk.vcf.gz $REGION \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/ld_continuous.py -v tissues/$TISSUE/$EGENE/vars.txt -f DS -a r \
    > tissues/$TISSUE/$EGENE/vars.cts.ld


