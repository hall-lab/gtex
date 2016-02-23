#!/bin/bash -e

if [[ -z $1 ]] || [[ -z $2 ]]
then
    echo "usage: $0 <LIST> <VCF>"
    exit 1
fi

LIST=$1
VCF=$2

# job index
EGENE=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 2`
TISSUE=`cat $LIST | sed -n ${LSB_JOBINDEX},${LSB_JOBINDEX}p | cut -f 3`

# echo -e "$LSB_JOBINDEX\t$EGENE\t$TISSUE"

mkdir -p tissues/$TISSUE/$EGENE

# ------------------------------------
# get expression
# EXPR=/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/rare_variants_2015-10-13/PEER.08052015/$TISSUE.peer.ztrans.txt.gz
EXPR=/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/depr/rare_variants_2015-08-07/peer_expr_z/$TISSUE.expr.peer.z.txt.gz
zcat $EXPR \
    | awk -v EGENE=$EGENE '{ if (NR==1) { for (i=2;i<=NF;++i) SAMPLES[i]=$i } if ($1==EGENE) { for (i=2;i<=NF;++i) print SAMPLES[i],SAMPLES[i],$i } }' OFS="\t" FS="\t" \
    | zjoin -w a -a stdin -b tissues/$TISSUE/keep.samples.txt -1 1 -2 1 \
    > tissues/$TISSUE/$EGENE/expr.phen

# ------------------------------------
# get SNV/indel variants

cat tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt \
    | awk -v PTHRES=$PTHRES -v EGENE=$EGENE '{ if ($1==EGENE && $2~"b37" && $4!~"nan") print $2,$4 }' OFS="\t" \
    | sort -k2,2g | head -n 1000 | cut -f 1 \
    > tissues/$TISSUE/$EGENE/top_1000_snv_indel.txt

# ------------------------------------
# get each sample's genotype and dosage for the SV

ESV=`cat tissues/$TISSUE/$EGENE/$EGENE.nom.eqtl.txt | awk '$4!~"nan$"' | awk '$2!~"_b37$"' | sort -k 4,4g | awk '{ print $2; exit }'`
REGION=`zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz | awk -v EGENE=$EGENE -v SLOP=1000000 '{ if ($1==EGENE) { POS_START=$3-SLOP; POS_END=$4+SLOP; if (POS_START<0) POS_START=0; print $2":"POS_START"-"POS_END } }'`

tabix -h $VCF $REGION \
    | vawk --header -v ESV=$ESV '{ if ($3==ESV) { if (I$SVTYPE=="CNV") { print $1,$2,$3,$4,$5,$6,$7,$8,$9,S$*$CN; print $1,$2,$3,$4,$5,$6,$7,$8,$9,S$*$DS; exit } else { print $1,$2,$3,$4,$5,$6,$7,$8,$9,S$*$GT; print $1,$2,$3,$4,$5,$6,$7,$8,$9,S$*$DS; exit } } }' \
    | grep -v '^##' \
    | cut -f 10- \
    | transpose \
    | awk -v ESV=$ESV '{ print $1,ESV,$2,$3 }' OFS='\t' \
    > tissues/$TISSUE/$EGENE/sv.genotypes.txt

# --------------------------------------------
# get normalized SV genotypes for each individual
MEAN=`cat tissues/$TISSUE/$EGENE/sv.genotypes.txt | cut -f 4 | zstats | awk '$1=="arith." { print $3 }'`
STDEV=`cat tissues/$TISSUE/$EGENE/sv.genotypes.txt | cut -f 4 | zstats | awk '$1=="stdev:" { print $2 }'`
cat tissues/$TISSUE/$EGENE/sv.genotypes.txt | awk -v MEAN=$MEAN -v STDEV=$STDEV '{ if (STDEV=="") { print $1,$1,0 } else { Z=($4-MEAN)/STDEV; print $1,$1,Z } }' OFS="\t" \
    > tissues/$TISSUE/$EGENE/sv.z.genotypes.txt

# --------------------------------------------
# get other covariates
zjoin -a tissues/$TISSUE/$EGENE/sv.z.genotypes.txt -b <(less ../tissues/$TISSUE/$TISSUE.covariates.txt \
    | transpose \
    | awk '{ print $1,$(NF-1),$NF }' OFS="\t") -1 1 -2 1 \
    | cut -f 1-3,5- \
    > tissues/$TISSUE/$EGENE/qcovar.txt

# --------------------------------------------
# make genetic relatedness matrices (GRM)
rm -r tissues/$TISSUE/$EGENE/mott
mkdir -p tissues/$TISSUE/$EGENE/mott
tabix -h $VCF $REGION \
    | /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/src/grm.py \
    -s tissues/$TISSUE/keep.samples.txt \
    -v tissues/$TISSUE/$EGENE/top_1000_snv_indel.txt \
    -f DS \
    -a mott \
    -o tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel

rm tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel.reml-no-constrain*

# --------------------------------------------
# run gcta (no constrain)
gcta64 \
    --reml-no-constrain \
    --thread-num 1 \
    --grm-gz tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel \
    --pheno tissues/$TISSUE/$EGENE/expr.phen \
    --qcovar tissues/$TISSUE/$EGENE/sv.z.genotypes.txt \
    --reml-est-fix \
    --reml-pred-rand \
    --out tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel.reml-no-constrain

# --------------------------------------------
# run gcta with extra covariates (no constrain)
gcta64 \
    --reml-no-constrain \
    --thread-num 1 \
    --grm-gz tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel \
    --pheno tissues/$TISSUE/$EGENE/expr.phen \
    --qcovar tissues/$TISSUE/$EGENE/qcovar.txt \
    --reml-est-fix \
    --reml-pred-rand \
    --out tissues/$TISSUE/$EGENE/mott/top_1000_snv_indel.extra_cov.reml-no-constrain
