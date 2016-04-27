#!/usr/bin/env Rscript

# 2015-12-02

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.


# print usage
usage <- function() {
  cat(
    'usage: prog <assoc> <ld> <gene> <title> <gwas_snp> <output>
author: Colby Chiang (colbychiang@wustl.edu)
    assoc        table of p-value assoc with trait
    ld           ld with index variant
    gene         gene file
    title        title of plot
    gwas snp     id of gwas snp
    output
')
    }

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
ASSOC_FILE <- args[1]
LD_FILE <- args[2]
GENE_FILE <- args[3]
TITLE <- args[4]
GWAS_SNP <- args[5]
OUTPUT <- args[6]

# Check input args
if (is.na(ASSOC_FILE)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(LD_FILE)) {
  usage()
  quit(save='no', status=1)
}

# setwd('/Users/cchiang/research/hall/projects/gtex/expts/phase1_2015-04-06/merged_2015-06-26/fastqtl_2015-10-11/gwas_2015-11-23/locus_zoom_test/ENSG00000040487.8')

library('LSD')

## par(mfrow=c(2,1))

# read input data
p <- read.table(ASSOC_FILE, header=T, stringsAsFactors=F)
p$logp <- -log(p$pval,10)
ld <- read.table(LD_FILE, stringsAsFactors=F)
colnames(ld) <- c('id', 'index_var', 'r')
ld$r <- as.numeric(ld$r)
genes <- read.table(GENE_FILE, stringsAsFactors=F)


## print(ld[which(ld[,1]=="9_124422403_G_A_b37"),])
## head(ld)

p <- merge(p, ld, by = "id")
## print(p$r)

# square the r value
p$r2 <- as.numeric(p$r) ** 2
p$abs_r <- abs(as.numeric(p$r))

col.res <- 10
cols <-  colorpalette('matlablike2', col.res)
p$col <- cols[cut(p$abs_r, breaks=col.res, labels=F)]
p <- p[order(p$abs_r, decreasing=F),]


# make the index SV special
sv <- p$index_var[1]
sv.index <- which(p$id==sv)[1]
p[sv.index,]$col="#FF1493"

## print(p[,11:13])


## quit()
## print(p[sv.index,])


## print(sv)


# plot
## pdf(OUTPUT, height=7, width=7)
png(OUTPUT, height=7, width=7, units='in', res=150)

xrange <- abs(max(p$end) - min(p$start))
xlimits <- c(min(p$start) - 0.05 * xrange, max(p$end) + 0.15 * xrange)
ylimits <- c(-1,max(p$logp)*1.05)

plot(logp~start, data=p, pch=21, col='black', bg=p$col, cex=1.5, main=paste(strwrap(TITLE, 90), collapse='\n'), cex.main=0.6, axes=F, xlim=xlimits, ylim=ylimits, ylab='-log(p)', xlab=paste0('Chromosome ', p$chrom[1], ' position (kb)'))
box()
axis(2)
axis(1, at=axTicks(1), labels=axTicks(1)/1e3 )

egene.y <- -0.4
forward.y <- -0.6
reverse.y <- -0.8

# plot the GWAS snp special
GWAS_SNP_VECTOR <- unlist(strsplit(GWAS_SNP, split=","))
points(logp~start, data=p[p$id %in% GWAS_SNP_VECTOR,], pch=23, col='black', bg=p[p$id %in% GWAS_SNP_VECTOR,]$col, cex=2, lwd=2)
## points(logp~start, data=p[p$id==GWAS_SNP,], pch=23, col='red', bg=p[p$id==GWAS_SNP,]$col, cex=2.5, lwd=2)

# draw genes
for (i in 1:nrow(genes)) {
  gene.col <- 'black'

  if (genes[i,5]=="+") {
    x.start <- genes[i,2]
    x.end <- genes[i,3]
    y.val <- forward.y
  }
  else if (genes[i,5]=="-") {
    x.start <- genes[i,3]
    x.end <- genes[i,2]
    y.val <- reverse.y
  }
  if (genes[i,6]==1) {
    gene.col <- 'red'
    y.val <- egene.y
  }
  
  arrows(x0=x.start, y0=y.val, x1=x.end, y1=y.val, length=0.05, col=gene.col, lwd=1)
}

# box SV
sv <- p[grepl('b37$', p$id)==0,]
pad <- 0.15
sv.col <- 'black'
## y.val <- sv$logp + pad
## arrows(x0=sv$start, y0=y.val, x1=sv$end, y1=y.val, length=0)
rect(xleft=sv$start, ybottom=sv$logp - pad, xright=sv$end, ytop=sv$logp + pad, border=sv.col)
text(sv$start, sv$logp + 2*pad, sv$id, cex=0.6)

legend('topright', rev(paste0((0:(col.res-1))/col.res, '-', (1:col.res)/col.res)), fill=rev(cols), title='abs(R)', cex=0.8)

dev.off()

## # -------------------------------------------
## plotTriMatrix <- function(x) {
##   ## clear lower triangle
##   x[lower.tri(x)] <- NA

##   ## calculate diag
##   nr <- nrow(x)
##   nc <- ncol(x)
##   d <- sqrt(nr^2 + nc^2)
##   d2 <- 0.5 * d

##   ## empty plot area
##   plot(NA, type="n", xlim=c(0, d), ylim=c(0, d), xlab="", ylab="", asp=1)

##   ## plot matrix and rotate 45
##   rasterImage(as.raster(x),
##               xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
##               interpolate=FALSE, angle=45, col=heat.colors(12))
## }

## ld2 <- ld**2

## plotTriMatrix(ld2)

## ld2.tri[lower.tri(ld2.tri)] <- NA
## image(ld2)


