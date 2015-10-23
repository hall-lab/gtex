#!/bin/bash -e

TISSUE=$1
CHROM=$2

mkdir -p tissues/$TISSUE/scaled

tabix -h ../gtex.joint_sv_gatk.vcf.gz $CHROM \
    | /gscmnt/gc2719/halllab/users/cchiang/src/sandbox/column_select.py \
	-c <(zcat tissues/$TISSUE/$TISSUE.expr.bed.gz | cut -f 5- | head -n 1 | tr '\t' '\n') \
	-l 9 \
	-p "##" \
    | vawk --header '$7=="PASS" || $7=="LOW"' \
    | vawk --header '{ delete ALLELES; NS=0; MS=0; for (i=10;i<=NF;++i) { if (I$SVTYPE!="CNV") { if ($i!~"^0/0" && $i!~"^\\./\\.") NS+=1; if ($i!~"^1/1" && $i!~"^\\./\\.") MS+=1; } else { split($i,GT,":"); ALLELES[GT[2]]+=1; } } if (I$SVTYPE=="CNV") { CN10=0; for (A in ALLELES) { if (ALLELES[A]>=10) CN10 += 1; } } if ((I$SVTYPE=="CNV" && CN10>=2) || (NS>=10 && MS>=10)) print $0 }' \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/scale_gt.py \
    | bgzip -c > tissues/$TISSUE/scaled/$CHROM.vcf.gz
