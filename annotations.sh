# 2015-05-13

# annotation

pwd
# /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations

# --------------------------------
# Gencode
curl -OL ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz

# make BED file of gene transcripts
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="transcript"' \
    | awk '{ gsub("^chr","",$1); gsub("[\";]","",$10); print $1,$4-1,$5,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.genes.bed.gz

# make BED file of gene with strands
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="transcript"' \
    | awk '{ gsub("^chr","",$1); gsub("[\";]","",$10); print $1,$4-1,$5,$10,$7}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.genes.with_strand.bed.gz

# make bed file of exons
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="exon"' \
    | awk '{ gsub("^chr", "", $1); gsub("[\";]", "", $10); print $1,$4-1,$5,$10,$7 }' OFS="\t" \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.exons.bed.gz

# make BED file of <= 5kb of TSS
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="gene"' \
    | awk -v DIST=5000 '{ gsub("^chr","",$1); if ($7=="+") { TSS=$4; PSTART=TSS-DIST/2-1; PEND=TSS } else if ($7=="-") { TSS=$5; PSTART=TSS-1; PEND=TSS+DIST/2 } gsub("[\";]","",$10); print $1,PSTART,PEND,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.5kb_upstream_TSS.bed.gz

# make BED file of <= 5kb of 3' end
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="gene"' \
    | awk -v DIST=5000 '{ gsub("^chr","",$1); if ($7=="+") { threePrime=$5; PSTART=threePrime-DIST/2-1; PEND=threePrime } else if ($7=="-") { threePrime=$4; PSTART=threePrime-1; PEND=threePrime+DIST/2 } gsub("[\";]","",$10); print $1,PSTART,PEND,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.5kb_downstream_3prime.bed.gz

# make BED file of <= 10kb of TSS
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="gene"' \
    | awk -v DIST=10000 '{ gsub("^chr","",$1); if ($7=="+") { TSS=$4; PSTART=TSS-DIST/2-1; PEND=TSS } else if ($7=="-") { TSS=$5; PSTART=TSS-1; PEND=TSS+DIST/2 } gsub("[\";]","",$10); print $1,PSTART,PEND,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.10kb_upstream_TSS.bed.gz

# make BED file of <= 10kb of 3' end
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="gene"' \
    | awk -v DIST=10000 '{ gsub("^chr","",$1); if ($7=="+") { threePrime=$5; PSTART=threePrime-DIST/2-1; PEND=threePrime } else if ($7=="-") { threePrime=$4; PSTART=threePrime-1; PEND=threePrime+DIST/2 } gsub("[\";]","",$10); print $1,PSTART,PEND,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.10kb_downstream_3prime.bed.gz


# make BED file of introns (in transcript but subtracting the exons)
zcat /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.annotation.gtf.gz \
    | awk '$3=="transcript"' \
    | awk '{ gsub("^chr","",$1); gsub("[\";]","",$10); print $1,$4-1,$5,$10}' OFS="\t" \
    | awk '{ if ($2<0) $2=0; print }' OFS="\t" \
    | bedtools subtract -a stdin -b /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.exons.bed.gz \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/gencode.v19.introns.bed.gz


# ------------------------------------
# WARNING: these are based soley on histone marks (not chip or FAIRE-seq data) and are now abandoned (2015-11-10)
# chromHMM tracks from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/
for TRACK in \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHepg2HMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHuvecHMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhekHMM.bed.gz \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed.gz
do
    echo $TRACK
    TRACKBASE=`basename $TRACK`
    curl -s $TRACK \
        | gzip -cdfq \
        | sed 's/^chr//g' \
        | sort -k1,1V -k2,2n -k3,3n \
        | bgzip -c \
        > $TRACKBASE
done

# merge into a single track of enhancers
zcat wgEncodeBroadHmm*.bed.gz \
    | awk '$4~"Enhancer"' \
    | sort -k1,1V -k2,2n -k3,3n | bedtools merge -c 4 -o distinct \
    | bgzip -c \
    > chromHmm.9_tissues.bed.gz

# --------------------------------------------
# VISTA enhancers
# get data
curl -s "http://enhancer.lbl.gov/cgi-bin/imagedb3.pl?page_size=20000;show=1;search.result=yes;page=1;form=search;search.form=no;action=search;search.sequence=1" | grep "^>Human" | awk '{ print $2 }' FS="|" | sed -e 's/[:\-]/ /g' | sed 's/^chr//g' | awk '{ print $1,$2,$3 }' OFS="\t" | sort -k1,1V -k2,2n -k3,3n | gzip -c > /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/vista.enhancers.bed.gz

less /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations/vista.enhancers.bed.gz | wc -l
# 1739


# --------------------------------------------
# Roadmap epigenomics
pwd
# /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations

mkdir -p roadmap


# --------------------------------------------
# DNAse hypersensitivity sites
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz | gzip -cdfq | awk '{ gsub("^chr","",$1); print }' OFS="\t" | bgzip -c > wgEncodeRegDnaseClusteredV3.bed.gz

# divide these into different scores
for SCORE in 300 500 750 1000
do
    zcat wgEncodeRegDnaseClusteredV3.bed.gz | awk -v SCORE=$SCORE '$5>=SCORE' | bgzip -c > wgEncodeRegDnaseClusteredV3.score_$SCORE.bed.gz
    tabix wgEncodeRegDnaseClusteredV3.score_$SCORE.bed.gz
done

# -------------------------------------------------
# transcription factor binding sites

# schema:
# 
#   Database: hg19    Primary Table: wgEncodeRegTfbsClusteredV3    Row Count: 4,380,444   Data last updated: 2013-07-21
# Format description: BED5+ with two fields having variable number of experiment IDs and values (none zero-valued)
# field example SQL type    description
# bin   585 smallint(5) unsigned    Indexing field to speed chromosome range queries.
# chrom chr1    varchar(255)    Reference sequence chromosome or scaffold
# chromStart    10073   int(10) unsigned    Start position in chromosome
# chromEnd  10329   int(10) unsigned    End position in chromosome
# name  ZBTB33  varchar(255)    Name of item
# score 354 int(10) unsigned    Score from 0-1000
# expCount  2   int(10) unsigned    Number of experiment values
# expNums   204,246 longblob    Comma separated list of experiment numbers
# expScores 354,138 longblob    Comma separated list of experiment scores

curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz \
    | gzip -cdfq | awk '{ gsub("^chr","",$1); print }' OFS="\t" | bgzip -c > wgEncodeRegTfbsClusteredV3.bed.gz

# divide these into different scores
for SCORE in 300 500 750 1000
do
    zcat wgEncodeRegTfbsClusteredV3.bed.gz | awk -v SCORE=$SCORE '$5>=SCORE' | bgzip -c > wgEncodeRegTfbsClusteredV3.score_$SCORE.bed.gz
    tabix wgEncodeRegTfbsClusteredV3.score_$SCORE.bed.gz
done


# ---------------------------------------------
# Mobile elements
# 2015-06-05
# generate repeat elements track (http://genome.ucsc.edu/cgi-bin/hgTables)
# less than 200 millidiv
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
     | gzip -cdfq \
     | awk '{ gsub("^chr", "", $6); if ($3<200) print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
     | sort -k1,1V -k2,2n -k3,3n \
     | bgzip -c > repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz
tabix -p bed repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz

# greater than or equal to 200 millidiv
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
     | gzip -cdfq \
     | awk '{ gsub("^chr", "", $6); if ($3>=200) print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
     | sort -k1,1V -k2,2n -k3,3n \
     | bgzip -c > repeatMasker.ancient.gte200millidiv.b37.sorted.bed.gz
tabix -p bed repeatMasker.ancient.gte200millidiv.b37.sorted.bed.gz

# extract only the SINEs LINEs and SVAs from the
zcat repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz | awk '$4~"LINE" || $4~"SINE" || $4~"SVA"' \
    | bgzip -c > repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz

# ----------------------------------------------
# ENCODE (PSU, Yue Lab) predicted enhancers
# http://promoter.bx.psu.edu/ENCODE/download.html
wget http://promoter.bx.psu.edu/ENCODE/predicted_enhancer_human.tar.gz
tar -zxvf predicted_enhancer_human.tar.gz

wget http://promoter.bx.psu.edu/ENCODE/predicted_promoter_human.tar.gz
tar -zxvf predicted_promoter_human.tar.gz

# ----------------------------------------------
# ENCODE enhancers
# https://www.encodeproject.org/data/annotations/

# Distal H3K27ac annotations (cell type specific) 
wget https://www.encodeproject.org/files/ENCFF786PWS/@@download/ENCFF786PWS.bigBed








# http://genome.ucsc.edu/ENCODE/downloads.html

# ----------------------------------
# ENCODE TFBS
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz


# -------------------------------------
# ENCODE genome segmentations (chromHMM) (current version 2016-01-29)
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/
# TSS   Bright Red  Predicted promoter region including TSS
# PF    Light Red   Predicted promoter flanking region
# E Orange  Predicted enhancer
# WE    Yellow  Predicted weak enhancer or open chromatin cis regulatory element
# CTCF  Blue    CTCF enriched element
# T Dark Green  Predicted transcribed region
# R Gray    Predicted Repressed or Low Activity region

mkdir -p encode.segmentation
cd encode.segmentation
for TRACK in http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmGm12878.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmH1hesc.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHepg2.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHuvec.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmK562.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedGm12878.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedHelas3.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedHepg2.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedHuvec.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationCombinedK562.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayGm12878.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayH1hesc.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayHelas3.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayHepg2.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayHuvec.bed.gz \
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayK562.bed.gz
do
    echo $TRACK
    TRACKBASE=`basename $TRACK`
    curl -s $TRACK \
        | gzip -cdfq \
        | sed 's/^chr//g' \
        | sort -k1,1V -k2,2n -k3,3n \
        | bgzip -c \
        > $TRACKBASE
done

# union of elements
mkdir -p union
for ELEMENT in CTCF T R WE PF TSS E
do
    echo $ELEMENT
    zcat wgEncodeAwgSegmentationCombined*.bed.gz \
        | awk -v ELEMENT=$ELEMENT '$4==ELEMENT' \
        | sort -k1,1V -k2,2n -k3,3n | bedtools merge -c 4 -o distinct \
        | bgzip -c \
        > union/wgEncodeAwgSegmentationCombined.union.$ELEMENT.bed.gz
done

# intersection of elements
mkdir -p intersection
for ELEMENT in CTCF T R WE PF TSS E
do
    for CELL in Gm12878 H1hesc Helas3 Hepg2 Huvec K562
    do
        zcat wgEncodeAwgSegmentationCombined$CELL.bed.gz \
            | awk -v ELEMENT=$ELEMENT -v CELL=$CELL '{ if ($4==ELEMENT) print $1,$2,$3,$4,CELL }' OFS="\t"
    done | sort -k1,1V -k2,2n -k3,3n | bedtools merge -c 4,5,5 -o distinct,distinct,count_distinct \
        | bgzip -c \
        > intersection/wgEncodeAwgSegmentationCombined.intersect.$ELEMENT.bed.gz
done

# get those in at least 2 tissues
MIN_TISSUES=2
for ELEMENT in CTCF T R WE PF TSS E
do
    zcat intersection/wgEncodeAwgSegmentationCombined.intersect.$ELEMENT.bed.gz \
        | awk -v MIN_TISSUES=$MIN_TISSUES '$NF>=2' \
        | bgzip -c \
        > intersection/wgEncodeAwgSegmentationCombined.intersect.min_tissues_$MIN_TISSUES.$ELEMENT.bed.gz
done

# make a track that combines weak and strong enhancers
zcat union/wgEncodeAwgSegmentationCombined.union.WE.bed.gz union/wgEncodeAwgSegmentationCombined.union.E.bed.gz \
    | sort -k1,1V -k2,2n -k3,3n \
    | bedtools merge -c 4 -o distinct \
    | bgzip -c \
    > union/wgEncodeAwgSegmentationCombined.union.WE_plus_E.bed.gz

# ---------------------------------------------
# oreganno literature curated enhancers
curl -s http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/oreganno.txt.gz \
    | gzip -cdfq \
    | awk '{ gsub("^chr", "", $2); print }' OFS="\t" \
    | cut -f 2- \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > oreganno_b37.bed.gz
    
    | less


# -------------------------------------------
# funseq:
# * 1kg.phase1.snp.bed.gz   (bed format)
#   Content : all 1KG phaseI SNVs in bed format. 
#   Columns : chromosome , SNVs start position (0-based), SNVs end position, MAF (minor allele frequency)
#   Purpose : to filter out variants against 1KG SNVs based on allele frequencies. 
# * 1kg.phase1.snp.bed.gz.tbi
#   Index file of 1kg.phase1.snp.bed.gz
# * ENCODE.annotation.gz   (bed format)
#   Content : compiled annotation files from ENCODE, Gencode v7 and others, includes DHS, TF peak, Pseudogene, ncRNA, enhancers
#   Columns : chromosome , annotation start position (0-based), annotation end position, annotation name.
#   Purpose : to find SNVs in annotated regions. 
# * ENCODE.tf.bound.union.bed  (bed format)
#   Content : transcription factor (TF) binding motifs under ENCODE TF peaks. 
#   Columns : chromosome, start position (0-based), end position, motif name, , strand, TF name
# * gene.strong.selection
#   Content : genes under strong negative selection (fraction of rare SNVs among non-synonymous variants). 
# * drm.gene.bed
#   Content : distal regulatory modules with gene information, generated with new algorithm (~769K elements with ~17K genes). 
#   Purpose : to associate noncoding SNVs with genes
# * conserved.bed
#   Content : regions defined as 'ultra-conserved' regions (Bejerano, et al., 2004). 
# * sensitive.bed
#   Content : sensitive regions defined in (Khurana, et al., 2013)
#   User can generate novel sensitive regions using scripts under '1. Building Data Context'. 
# * All_hg19_RS.bw
#   File downloaded from UCSC. Gerp score file. 
# * hot.regions.bed
#   Content : highly occupied regions defined in (Yip, et al., 2012)
# * motif.PFM
#   Content : position frequency matrix
#   Purpose : used for motif breaking & motif gaining calculation
# * motif.score.cut 
#   File used to speed up the motif-gaining analysis.
#   Can be generated by '5.PWM.score.cut.pl' under '1. Building Data Context'. 
# * regulatory.network
#   Content : two columns, (TF, genes_regulated_by_TF)
#   Purpose : find TFs regulating known cancer genes. 
# * human_g1k_v37.fasta
#   Human reference genome (hg19)
# * human_g1k_v37.fasta.fai
#   Indexed reference genome (hg19)
# * human_ancestor_GRCh37_e59.fa
#   Human ancestral genome (hg19) 
#   Purpose : for motif breaking calculation in personal or germ-line genome. 
#   * Note :  for only somatic analysis, these files are not needed. 
# * human_ancestor_GRCh37_e59.fa.fai
#   Indexed ancestral genome (hg19) 
#   * Note :  for only somatic analysis, these files are not needed.
# * weighted.score.txt
#   Weighted scores for features

# sensitive regions
curl -s http://archive.gersteinlab.org/funseq2.1.0_data/sensitive.bed \
    | awk '{ gsub("^chr", "", $1); print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > funseq-2.1.0.sensitive.bed.gz

# TF bound
curl -s http://archive.gersteinlab.org/funseq2.1.0_data/ENCODE.tf.bound.union.bed \
    | awk '{ gsub("^chr", "", $1); print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > funseq-2.1.0.ENCODE.tf.bound.union.bed.gz

# ultra conserved
curl -s http://archive.gersteinlab.org/funseq2.1.0_data/conserved.bed \
    | awk '{ gsub("^chr", "", $1); print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > funseq-2.1.0.conserved.bed.gz

# hot regions
curl -s http://archive.gersteinlab.org/funseq2.1.0_data/hot.regions.bed \
    | awk '{ gsub("^chr", "", $1); print }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > funseq-2.1.0.hot.regions.bed.gz


# -----------------------------------------------------
# Dragon Enhancers Database

pwd
# /gscmnt/gc2719/halllab/users/cchiang/projects/gtex/annotations

mkdir -p DENdb
cd DENdb

wget http://www.cbrc.kaust.edu.sa/dendb/src/ChIP_seq_TF.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/dnase.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/enhancer_dnase.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/enhancers.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/enhancer_interactions.tsv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/enhancers_targets.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/FANTOM_expression.csv.zip
wget http://www.cbrc.kaust.edu.sa/dendb/src/tf.csv.zip

# unzip
ls *.zip | xargs -I{} unzip {}

for FILE in `ls *.csv`
do
    FBASE=`basename $FILE .csv`
    cat $FILE \
        | tr ',' '\t' \
        | awk '{ gsub("^chr","",$2); print $2,$3,$4,$1,$5,$6 }' OFS="\t" \
        | sort -k1,1V -k2,2n -k3,3n \
        | bedtools merge -c 6,6 -o distinct,count_distinct \
        | bgzip -c \
        > $FBASE.bed.gz
done

# slice the enhancers
for SCORE in {2..5}
do
    zcat enhancers.bed.gz | awk -v SCORE=$SCORE '$5>=SCORE' | bgzip -c > enhancers.score_$SCORE.bed.gz
done

# slice the dnase
for SCORE in {2..5}
do
    zcat dnase.bed.gz | awk -v SCORE=$SCORE '$5>=SCORE' | bgzip -c > dnase.score_$SCORE.bed.gz
done

# --------------------------------------------------
# 2016-01-22

# STR track, since our BNDs seem to be enriched in STRs
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
    | gzip -cdfq \
    | awk '{ gsub("^chr", "", $6); print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
    | awk '$4~"^Simple_repeat"' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > repeatMasker.simple_repeat.b37.sorted.bed.gz
tabix -p bed repeatMasker.simple_repeat.b37.sorted.bed.gz

# STR plus LTR track
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
    | gzip -cdfq \
    | awk '{ gsub("^chr", "", $6); print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
    | awk '$4~"^Simple_repeat" || $4~"^LTR"' \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > repeatMasker.simple_repeat_ltr.b37.sorted.bed.gz
tabix -p bed repeatMasker.simple_repeat_ltr.b37.sorted.bed.gz

# ---------------------
# better STR track
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
zcat simpleRepeat.txt.gz \
    | awk '{ gsub("^chr", "", $2); print }' OFS="\t" \
    | cut -f 2- \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > simpleRepeat.bed.gz
tabix -p bed simpleRepeat.bed.gz
zcat simpleRepeat.bed.gz | bedtools merge | awk '{ print $3-$2 }' | zsum
# 73061420


# ------------------------------------------------
# seg dups

curl -s http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz | gzip -cdfq \
    | awk '{ gsub("^chr","",$2); gsub("^chr","",$5) ; print $2,$3,$4,$5,$6,$7 }' OFS="\t" \
    | sort -k1,1V -k2,2n -k3,3n \
    | bgzip -c \
    > segdups.bed.gz


