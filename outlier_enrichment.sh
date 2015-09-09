#!/bin/bash -e

USAGE='\tusage: '$0' <VCF> <outliers> <distance> <cadd_min>'

if [ $# -lt 3 ]
then
    echo -e "\n$USAGE\n"
    exit
fi

CADD_MIN=0

VCF=$1
OUTLIERS=$2
DISTANCE=$3

if [ ! -z $4 ]
then
    CADD_MIN=$4
fi

# only look for private and doubleton
NSAMP_MAX=2

while read line
do
    # SAMPLE=`echo "$line" | cut -f 9`
    BED_REGION=`echo "$line" | vawk -v DISTANCE=$DISTANCE '{ if (I$SVTYPE~"^DEL" || I$SVTYPE=="DUP") BED=$1"\t"$2-DISTANCE"\t"I$END+DISTANCE; else if (I$SVTYPE=="INV") { BED=$1"\t"$2-DISTANCE"\t"$2+1+DISTANCE; BED=$1"\t"I$END-1-DISTANCE"\t"I$END+DISTANCE } else if (I$SVTYPE=="BND") BED=$1"\t"$2-1-DISTANCE"\t"$2+DISTANCE; print BED,$3,I$SVTYPE,$9 }'`
    echo -e "$BED_REGION" \
	| bedtools intersect -a /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.genes.bed.gz -b stdin -wo \
	| zapdups -u -k 4,8,10
done < <(zless $VCF | vawk --header -v NSAMP_MAX=$NSAMP_MAX -v CADD_MIN=$CADD_MIN 'I$NSAMP<=NSAMP_MAX && I$SVSCORE_NOTRUNC>CADD_MIN' | awk '{ if ($0~"^##") next; if ($0~"^#CHROM") { for (i=10;i<=NF;++i) COLSAMP[i]=$i } else { for (i=10;i<=NF;++i) { if ($i~"0/1" || $i~"1/1") print $1,$2,$3,$4,$5,$6,$7,$8,COLSAMP[i] } } }' OFS="\t")
