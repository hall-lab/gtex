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



# --------------------------------------------
# DNAse hypersensitivity sites
curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz | gzip -cdfq | awk '{ gsub("^chr","",$1); print }' OFS="\t" | bgzip -c > wgEncodeRegDnaseClusteredV3.bed.gz
curl -OL http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredSources.tab

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
# ENCODE genome segmentations (chromHMM)
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

curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz \
     | gzip -cdfq \
     | awk '{ gsub("^chr", "", $6); if ($3<200) print $6,$7,$8,$12"|"$13"|"$11,$3,$10 }' OFS="\t" \
     | sort -k1,1V -k2,2n -k3,3n \
     | bgzip -c > repeatMasker.recent.lt200millidiv.b37.sorted.bed.gz









