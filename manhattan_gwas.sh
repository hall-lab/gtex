#!/bin/bash -e

if [[ -z $1 ]] || [[ -z $2 ]]
then
    echo "usage: $0 <LIST> <OUTDIR> <MINLD>"
    exit 1
fi

LIST=$1
OUTDIR=$2
MINLD=$3

if [[ -z "$MINLD" ]]
then
    MINLD=0
fi

# LSB_JOBINDEX=20

# job index
TITLE=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 1`
EGENE=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 2`
TISSUE=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 3`
ASSAY=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 4`
GWAS_SNP=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 5`

# OUTDIR=plots
mkdir -p $OUTDIR

echo -e "$LSB_JOBINDEX\t$EGENE\t$TISSUE"

NOM_ASSOC_FILE=../caviar/tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt
# echo $NOM_ASSOC_FILE

# ------------------------------------
# make directory

 mkdir -p tissues/$TISSUE/$EGENE

# ---------------------------------
# get variants
cat $NOM_ASSOC_FILE \
    | awk '{ if ($4!~"nan") print $2,$4 }' OFS="\t" \
    | sort -k2,2g | cut -f 1 \
    > tissues/$TISSUE/$EGENE/vars.txt

# -----------------------------------
# get variant locations
REGION=`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz | awk -v EGENE=$EGENE -v SLOP=10000000 '{ if ($1==EGENE) { POS_START=$3-SLOP; POS_END=$4+SLOP; if (POS_START<0) POS_START=0; print $2":"POS_START"-"POS_END } }'`

tabix ../tissues/$TISSUE/$TISSUE.joint_sv_gatk.scaled.sup_10_samp.high_conf.vcf.gz $REGION \
    | zjoin -a stdin -b tissues/$TISSUE/$EGENE/vars.txt -1 3 -2 1 \
    | vcfToBed | cut -f -6 \
    > tissues/$TISSUE/$EGENE/var_locs.bed

# ----------------------------------
# get association with the trait
zjoin -a tissues/$TISSUE/$EGENE/var_locs.bed -b $NOM_ASSOC_FILE -1 4 -2 2 \
    | awk 'BEGIN { print "chrom\tstart\tend\tid\tref\talt\tbeta\tpval" } { print $1,$2,$3,$4,$5,$6,$12,$11 }' OFS="\t" \
    | zapdups -u -k 4 \
    > tissues/$TISSUE/$EGENE/var_assoc.bed

# -------------------------------
# get LD between the SV and the other variants
SV_ID=`cat $NOM_ASSOC_FILE | awk '$2!~"_b37$"' | awk '{ if ($4!~"nan") print $2,$5 }' OFS="\t" | sort -k2,2g | cut -f 1 | head -n 1`
echo $SV_ID
# REGION=`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz | awk -v EGENE=$EGENE -v SLOP=10000000 '{ if ($1==EGENE) { POS_START=$3-SLOP; POS_END=$4+SLOP; if (POS_START<0) POS_START=0; print $2":"POS_START"-"POS_END } }'`

tabix -h /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/merged_2015-06-26/gtex.joint_sv_gatk.vcf.gz $REGION \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/ld_continuous.py -v tissues/$TISSUE/$EGENE/vars.txt -f DS -a r -I $SV_ID \
    > tissues/$TISSUE/$EGENE/vars.cts.ld

# --------------------------------------
# get the genes

tabix /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.genes.with_strand.bed.gz $REGION \
    | awk '{ print $4"_"$5"_"$1,$2,$3,$4 }' OFS="\t" | sort -k1,1V -k2,2n -k3,3 \
    | bedtools merge \
    | awk -v EGENE=$EGENE '{ IS_EGENE=0; split($1,x,"_"); if (x[1]==EGENE) IS_EGENE=1; print x[3],$2,$3,x[1],x[2],IS_EGENE }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3 \
    > tissues/$TISSUE/$EGENE/genes.txt

# # -----------------------------------

GWAS_LD=`cat tissues/$TISSUE/$EGENE/vars.cts.ld | awk -v GWAS_SNP=$GWAS_SNP '$1==GWAS_SNP { print $3 }'`
PASS=`awk -v MINLD=$MINLD -v GWAS_LD=$GWAS_LD 'BEGIN { if (GWAS_LD*GWAS_LD >= MINLD*MINLD) print 1; else print 0; exit }'`

# # echo $GWAS_LD2 $MINLD2
# echo $MINLD $GWAS_LD
# echo $PASS

if [[ "$PASS" == 0 ]]
then
    exit 0
fi

# ----------------------------------
# draw the plot

# GWAS_SNP=`cat tissues/$TISSUE/$EGENE/gwas_var.txt | cut -f 1`
RANK=`printf "%04d\n" $LSB_JOBINDEX`

echo $TISSUE $EGENE
/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/manhattan_plot.R \
    tissues/$TISSUE/$EGENE/var_assoc.bed \
    tissues/$TISSUE/$EGENE/vars.cts.ld \
    tissues/$TISSUE/$EGENE/genes.txt \
    "$RANK: $TITLE" \
    $GWAS_SNP \
    $OUTDIR/$RANK.$TISSUE.$EGENE.locus.png
