#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
out <- args[2]

d = read.table(file, header=F, stringsAsFactors=F)

## colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "ppval", "bpval")
colnames(d) = c("pid", "sid", "dist", "beta", "npval", "ppval", "bpval")

png(out, width=8, height=8, res=300, units='in')
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot", ylim=c(0,1))
abline(0, 1, col="red")
dev.off()

