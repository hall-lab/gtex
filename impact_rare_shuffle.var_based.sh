#!/usr/bin/env bash

PTILE=$1
SHUFFLED_DIR=$2

SVSCORE_FILE=/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/merged_2015-06-26/data/svscore.notrunc.txt.gz

for i in {1..1000}
do
    SHUFFLED_FILE=$SHUFFLED_DIR/shuffle.$i.bed.gz

    NCOL=`zcat $SHUFFLED_FILE | awk '{ print NF; exit }'`
    zjoin -r a -a <(zcat $SHUFFLED_FILE | awk '{ gsub("_[12]$","",$4); print }' OFS="\t") -b <(zcat $SVSCORE_FILE | awk '{ gsub("_[12]$","",$1); print }' OFS="\t" | zapdups -u -k1 ) -1 4 -2 1 \
        | awk '{ if ($4!~"b37$") $9=$9";SVSCORE_NOTRUNC="$NF; print }' OFS="\t" \
        | cut -f -$NCOL \
        | vawk -c 9 -v RUN=$i '{ split($4,id,"_"); if ($4~"b37$") { if (length(id[3])==1 && length(id[4])==1) TYPE="SNV"; else TYPE="INDEL" } else TYPE="SV"; if (TYPE=="SV") SCORE=I$SVSCORE_NOTRUNC; else SCORE=I$CADD; print RUN,$4,TYPE,SCORE,$13,$14,$15,$16,$17 }' \
        | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/impact_score_cdf.py --sv data/sv.cadd.percentile.txt  --snv data/snv.cadd.percentile.txt --indel data/indel.cadd.percentile.txt \
        | sort -k10,10gr \
        | zapdups -u -k 2,3 \
        | awk -v RUN=$i -v PTILE=$PTILE 'BEGIN { SV=0; SNV=0; INDEL=0; } { if ($10<PTILE) next; if ($3=="SV") SV+=1; else if ($3=="SNV") SNV+=1; else if ($3=="INDEL") INDEL+=1 } END { print PTILE,RUN,SV,SNV,INDEL }' OFS="\t"
done
