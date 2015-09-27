#!/usr/bin/env Rscript

# 2015-09-24

# print usage
usage <- function() {
    cat(
          'usage: expression.heatmap.R <tissue>
expression.heatmap.R
author: Colby Chiang (colbychiang@wustl.edu)
    matrix       expression matrix
')
  }

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
FILE <- args[1]

# Check input args
if (is.na(FILE)) {
    usage()
      quit(save='no', status=1)
  }


library(gplots)

# number of variants to downsample to
num.ds <- 1000

# simulate the colors of Matlab
myCols=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

# set your display range here
pairs.breaks <- seq(-3, 3, by=.05)
## pairs.breaks
## length(pairs.breaks)

mat <- read.table(FILE, stringsAsFactors=FALSE, header=FALSE)

colnames(mat) <- mat[1,]
rownames(mat) <- mat[,1]
mat <- as.matrix(mat[-1,-1])
class(mat) <- 'numeric'
mat.ds <- mat[seq(1, nrow(mat), ceiling(nrow(mat)/num.ds)),]

## prior to PEER correction
png(paste0(FILE, '.expr.clust.png'), width=8, height=8, res=150, units='in')
heatmap.2(as.matrix(mat.ds), distfun=function(x) dist(x,method='euclidian'), hclustfun=function(x) hclust(x,method='ward.D2'),
                              trace='none', dendrogram='both', Rowv=TRUE, Colv=TRUE, breaks=pairs.breaks, col=colorRampPalette(myCols), scale='none', symkey=T, na.color='grey', density.info='histogram', cexRow=0.2, cexCol=0.5, main=paste0(FILE, '\nexpr clustering'))
dev.off()
