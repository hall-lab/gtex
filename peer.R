#!/usr/bin/env Rscript

# 2015-07-15

# print usage
usage <- function() {
  cat(
    'usage: peer.R <tissue>
peer.R
author: Colby Chiang (cc2qe@virginia.edu)
  tissue             tissue name
')
}

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
TISSUE <- args[1]

# Check input args
if (is.na(TISSUE)) {
  usage()
  quit(save='no', status=1)
}

# expression z score heat maps
# data here from David Deluca's 2015-01-12 eQTL input (log normalized) values

# setwd('/Users/cchiang/research/hall/projects/gtex/expts/phase1_2015-04-06/rare_variants_2015-07-14/R')

# -------------------------------------------
# SETUP

library(gplots)
library(peer)

# number of variants to downsample to
num.ds <- 1000

# simulate the colors of Matlab
myCols=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

# set your display range here
pairs.breaks <- seq(-3, 3, by=.05)
## pairs.breaks
## length(pairs.breaks)

OUT.PREFIX <- 'peer_expr'

# ----------------------------------

## TISSUE <- 'Heart_Left_Ventricle'
## TISSUE <- 'Whole_Blood'
## TISSUE <- 'test'
## TISSUE <- 'Brain_Hippocampus'

EXPR.DIR <- '/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_geneLevelNormalizedExpressionMatrices'

mat <- read.table(paste0(EXPR.DIR, '/', TISSUE, '_Analysis.expr.txt.gz'), stringsAsFactors=FALSE, header=FALSE)

colnames(mat) <- mat[1,]
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[-1,-1])
class(mat) <- 'numeric'
mat.ds <- mat[seq(1, nrow(mat), ceiling(nrow(mat)/num.ds)),]

## prior to PEER correction
## png(paste0(TISSUE, '.expr.clust.png'), width=8, height=8, res=150, units='in')
## heatmap.2(as.matrix(mat.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
          ## trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
## dev.off()


## ----------------------------------------
# peer correct

## dim(mat)

model <- PEER()
PEER_setPhenoMean(model, t(as.matrix(mat)))

dim(PEER_getPhenoMean(model))

# set covariates
COV.DIR <- '/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_covariates'
covs = read.table(paste0(COV.DIR, '/', TISSUE, '_Analysis.covariates.txt'), header=TRUE)
rownames(covs) <- covs[,1]
covs <- covs[,-1]
# only use gender
covs <- covs[which(rownames(covs)=="gender"),]
PEER_setCovariates(model, t(as.matrix(covs)))

# set number of peer factors
## N < 150, use 15  PEERs, 150<=N<250, use 30 PEERs, N >=250 use 35 PEERs
if (ncol(mat) < 150) {
  numcov <- 15
} else if (ncol(mat) < 250) {
  numcov <- 30
} else if (ncol(mat) >= 250) {
  numcov <- 35
}

PEER_setNk(model, numcov)
PEER_getNk(model)

## PEER_setTolerance(model, 0.1)
## PEER_setVarTolerance(model, 0.001)

PEER_update(model)

# diag
pdf(paste0(OUT.PREFIX, '/', TISSUE, '.peer.diag.pdf'), width=6, height=8)
PEER_plotModel(model)
dev.off()
## cor(PEER_getCovariates(model)[,1], PEER_getX(model)[,2])


factors = t(PEER_getX(model))
weights = PEER_getW(model)
precision = PEER_getAlpha(model)

residuals = t(PEER_getResiduals(model))
rownames(residuals) <- rownames(mat)
colnames(residuals) <- colnames(mat)

residuals.ds <- residuals[seq(1, nrow(residuals), ceiling(nrow(residuals)/num.ds)),]

png(paste0(OUT.PREFIX, '/', TISSUE, '.expr.peer.clust.png'), width=8, height=8, res=150, units='in')
heatmap.2(as.matrix(residuals.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
          trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(TISSUE, '\nexpr clustering'))
dev.off()

gz1 <- gzfile(paste0(OUT.PREFIX, '/', TISSUE, ".expr.peer.txt.gz"), "w")
write.table(cbind(rownames(residuals), residuals), file=gz1, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
close(gz1)
