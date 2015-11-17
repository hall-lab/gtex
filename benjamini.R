#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
out <- args[2]
fdr <- args[3]

d = read.table(file, header=F, stringsAsFactors=F)


## colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "ppval", "bpval")
colnames(d) = c("pid", "sid", "dist", "beta", "npval", "ppval", "bpval")

d$bh = p.adjust(d$bpval, method="fdr")
write.table(d[which(d$bh <= fdr & d$bpval>0),], out, quote=F, row.names=F, col.names=T, sep="\t")

