#!/bin/bash -e

TISSUE=$1
CHROM=$2

tabix -h ../gtex.joint_sv_gatk.vcf.gz $CHROM \
    | /gscmnt/gc2719/halllab/users/cchiang/src/sandbox/column_select.py \
	-c <(zcat tissues/$TISSUE/$TISSUE.expr.bed.gz | cut -f 5- | head -n 1 | tr '\t' '\n') \
	-l 9 \
	-p "##" \
    | vawk --header '$7=="PASS"' \
    | vawk --header 'I$AF>=0.05 && I$AF<=0.95' \
    | /gscmnt/gc2719/halllab/users/cchiang/src/svtyper/scripts/vcf_modify_header.py -i DS -c FORMAT -t Float -n A -d "Dosage of alt allele" \
    | awk '{ if ($0~"^#" || ($3!~"^LUMPY" && $3!~"^GS")) { print; next } $9=$9":DS"; split($9,fmt,":"); if ($3~"^LUMPY") { field="AB"; scalar=2 } if ($3~"^GS_") { field="CN"; scalar=1 } for (i=1;i<=length(fmt);++i) { if (fmt[i]==field) fmt_idx=i } for (i=10;i<=NF;++i) { split($i,gt,":"); $i=$i":"gt[fmt_idx]*scalar } print }' OFS="\t" \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/scale_lumpy_gt.py


