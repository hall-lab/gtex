#!/usr/bin/env Rscript

# ------------------
# functions
resid <- function(x, y) {
  res <- abs(y-x)/sqrt(2)
  return(res)
}
# ------------------

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
out <- args[2]
basename <- basename(file)

## d = read.table(file, header=F, stringsAsFactors=F)
d <- read.table(file, sep=' ', col.names=c('phenID', 'varID', 'distance', 'unknown', 'nomP', 'directP', 'betaP'), stringsAsFactors=F)
directP <- d$directP
betaP <- d$betaP

## Use densCols() output to get density at each point
x <- densCols(directP, betaP, colramp=colorRampPalette(c("black", "white")))
d$dens <- col2rgb(x)[1,] + 1L

## Map densities to colors
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100", "#9B0000"))(256)
d$col <- cols[d$dens]

## Plot it, reordering rows so that densest points are plotted on top
# pdf(paste0(basename(args[1]), '.pdf'), height=5, width=5)
png(out, height=5, width=5, units='in', res=300)
par(mar=c(4,4,3,1))
plot(betaP~directP, data=d[order(d$dens),], pch=20, col=col, cex=1.5, xlab='Direct method', ylab='Beta approximation', xlim=c(0,1), ylim=c(0,1), main=basename)

max.resid <- 0.1
abline(a=0, b=1, col='red')
abline(a=-max.resid * sqrt(2), b=1)
abline(a=max.resid * sqrt(2), b=1)

num.outliers <- 0
r <- resid(d$directP, d$betaP)
num.outliers <- length(r[which(r > max.resid)])

inlier.pct <- 1 - (num.outliers / nrow(d))
text(0.2, 0.8, paste0('Between black lines:\n', round(inlier.pct,6)))
dev.off()

