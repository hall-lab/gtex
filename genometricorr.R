#!/usr/bin/env Rscript

# print usage
usage <- function() {
    cat(
    'usage: rcorr.R <file>
    rcorr.R
    author: Colby Chiang (cc2qe@virginia.edu)
    description: calculate R and R^2 values from two columns of numbers
        outputs two tab-delimited columns: R and R^2
    positional arguments:
        file        File with two columns of numerical values,
                    tab-delimited.
    ')
}


# source("http://bioconductor.org/biocLite.R")
# biocLite("IRanges")
# biocLite("GenomicRanges")
# biocLite("GenomicFeatures")
# biocLite("rtracklayer")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# install.packages("gplots")
# install.packages("RCurl")

# install.packages('GenometriCorr',repos='http://genometricorr.sourceforge.net/R/',type='source')

library('GenometriCorr')

# test <- import('/Users/cchiang/research/hall/projects/gtex/expts/phase1_2015-04-06/merged_2015-06-26/fastqtl_2016-01-04_low/feature_enrichment_2016-01-06/cav1.bed.gz')
# feature <- import('/Users/cchiang/research/hall/projects/gtex/annotations/encode.segmentation/union/wgEncodeAwgSegmentationCombined.union.E.bed.gz')

test.file <- args[1]
feature.file <- args[2]
test <- import(test.file)
feature <- import(feature.file)

# VisualiseTwoIRanges(ranges(test[seqnames(test) == "1"]), ranges(feature[seqnames(feature) == "1"]), nameA = "test", nameB = "feature", chrom_length = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)["chr1"], title = "test and feature on chr1 of Hg19")

pn.area <- 100
pn.dist <- 100
pn.jacc <- 100
gc <- GenometriCorrelation(test, feature, ecdf.area.permut.number=pn.area, mean.distance.permut.number=pn.dist, jaccard.measure.permut.number=pn.jacc, keep.distributions=TRUE, showProgressBar=FALSE)


print(gc)
graphical.report(gc, pdffile='~/Desktop/test.pdf')
visualize(gc, pdffile='~/Desktop/test2.pdf')



read.table('/Users/cchiang/research/genomes/GRCh37/human_g1k_v37.sizes.txt')