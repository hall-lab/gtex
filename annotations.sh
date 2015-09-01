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
    curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz \
        | gzip -cdfq \
        | sed 's/^chr//g' \
        | sort -k1,1V -k2,2n -k3,3n \
        | bgzip -c \
        > $TRACKBASE
done

# merge into a single track of enhancers
zcat chromHMM_tissues/wgEncodeBroadHmm*.bed.gz \
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












