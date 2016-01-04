#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
out <- args[2]
fdr <- as.numeric(args[3])

d = read.table(file, header=F, stringsAsFactors=F)

colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "beta", "ppval", "bpval")
## colnames(d) = c("pid", "sid", "dist", "beta", "npval", "ppval", "bpval")

d$bonferroni = p.adjust(d$bpval, method="bonferroni")
write.table(d[which(d$bonferroni <= fdr & d$bpval>0),], out, quote=F, row.names=F, col.names=T, sep="\t")

