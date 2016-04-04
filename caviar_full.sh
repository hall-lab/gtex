#!/bin/bash -e

if [[ -z $1 ]]
then
    echo "usage: $0 <LIST>"
    exit 1
fi

LIST=$1

# job index
EGENE=`zless $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 1`
TISSUE=`zless $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 2`

echo -e "$LSB_JOBINDEX\t$EGENE\t$TISSUE"

# TISSUE=Whole_Blood
# EGENE=ENSG00000271523.1

# ---------------------------------------
# get variants
less tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt \
    | awk -v EGENE=$EGENE '{ EVENT=$2; if ($2~"LUMPY_BND") { gsub("_[12]$","",EVENT);} if ($1==EGENE && $2~"_b37$" && $5!~"nan") print $0,EVENT }' OFS="\t" \
    | sort -k5,5g | zapdups -u -k 7 | head -n 100 | cut -f 2 \
    > tissues/$TISSUE/$EGENE/top_100_snv_indel.txt
less tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt \
    | awk -v EGENE=$EGENE '{ EVENT=$2; if ($2~"LUMPY_BND") { gsub("_[12]$","",EVENT);} if ($1==EGENE && $2!~"_b37$" && $5!~"nan") print $0,EVENT }' OFS="\t" \
    | sort -k5,5g | zapdups -u -k 7 | head -n 1 | cut -f 2 \
    | cat - tissues/$TISSUE/$EGENE/top_100_snv_indel.txt \
    > tissues/$TISSUE/$EGENE/vars.txt

# ---------------------------------------
# get t-stat
cat tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt \
    | tr '\t' ' ' \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/fastqtl_qvalue.py \
        -p ../tissues/$TISSUE/$TISSUE.joint_sv_gatk.scaled.eqtl.txt.gz \
        -f matrixeqtl \
    | sed 1d \
    | zjoin -a stdin -b tissues/$TISSUE/$EGENE/vars.txt -1 1 -2 1 \
    | cut -f 1,4 \
    > tissues/$TISSUE/$EGENE/t-stat.txt

# ---------------------------------------
# get LD
REGION=`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz | awk -v EGENE=$EGENE -v SLOP=10000000 '{ if (\$1==EGENE) { POS_START=\$3-SLOP; POS_END=\$4+SLOP; if (POS_START<0) POS_START=0; print \$2":"POS_START"-"POS_END } }'`;

tabix -h ../../gtex.joint_sv_gatk.vcf.gz $REGION \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/ld_continuous.py -v tissues/$TISSUE/$EGENE/vars.txt -f DS -a r \
    > tissues/$TISSUE/$EGENE/vars.cts.ld

# ---------------------------------------
# run caviar
CAVIAR -o tissues/$TISSUE/$EGENE/caviar.t-stat.cts -l tissues/$TISSUE/$EGENE/vars.cts.ld -z tissues/$TISSUE/$EGENE/t-stat.txt -c 1


