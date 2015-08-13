#!/usr/bin/env Rscript

# 2015-08-11

# mahalanobis distance test

# setwd('/Users/cchiang/research/hall/projects/gtex/expts/phase1_2015-04-06/rare_variants_2015-08-07/mahalanobis_dist_example_2015-08-11/R')

# set seed for reproducibility
set.seed(10)


suppressMessages(library('gplots'))

imp.col <- function (a) {
	imputed <- a
	for (j in 1:ncol(a)) {
		missing <- is.na(a[,j])
		n.missing <- sum(missing)
		obs <- a[,j][!missing]
		imputed[,j][missing] <- sample(obs, n.missing, replace=TRUE)
	}
	return(imputed)
}

imp.colmean <- function (a) {
	imputed <- a
	for (j in 1:ncol(a)) {
		missing <- is.na(a[,j])
		n.missing <- sum(missing)
		obs <- a[,j][!missing]
		imputed[,j][missing] <- mean(obs)
	}
	return(imputed)
}

imp.row <- function (a) {
	imputed <- a
	for (i in 1:nrow(a)) {
		missing <- is.na(a[i,])
		n.missing <- sum(missing)
		obs <- a[i,][!missing]
		imputed[i,][missing] <- sample(obs, n.missing, replace=TRUE)
	}
	return(imputed)
}

imp.rowmean <- function (a) {
	imputed <- a
	for (i in 1:nrow(a)) {
		missing <- is.na(a[i,])
		n.missing <- sum(missing)
		obs <- a[i,][!missing]
		imputed[i,][missing] <- mean(obs)
	}
	return(imputed)
}



min.nsamp <- 10
tissue.list <- c("Whole_Blood", "Cells_Transformed_fibroblasts", "Muscle_Skeletal", "Lung", "Artery_Tibial", "Adipose_Subcutaneous", "Thyroid", "Esophagus_Mucosa", "Skin_Sun_Exposed_Lower_leg", "Nerve_Tibial", "Esophagus_Muscularis", "Artery_Aorta", "Heart_Left_Ventricle")

f <- file('stdin')
open(f)
# full <- read.table(file, header=T, row.names=NULL, stringsAsFactors=F)
prev.gene <- NULL
x <- NULL


# read header
header <- unlist(strsplit(readLines(f,n=1), '\t'))

# for (i in 1:nrow(full)) {
while(length(line <- readLines(f,n=1)) > 0) {
	v <- matrix(unlist(strsplit(line, '\t')), nrow=1)
	gene <- v[1]

	if (!is.null(prev.gene) && gene != prev.gene) {
		rownames(x) <- x[,2]
		colnames(x) <- header
		
		x <- x[,-(1:2)]
		x <- as.matrix(x)
		suppressWarnings(class(x) <- 'numeric')
		x <- t(x)
		# x <- x[which(apply(x, 1, function(z) sum(!is.na(z))) >= min.nsamp),]
		x <- x[, which(colnames(x) %in% tissue.list)]
                print(x)
                print(which(colnames(x) %in% tissue.list))
                print(ncol(x))
                print(is.null(x))
                if (!is.null(x) && ncol(x) > 0) {
                    x <- imp.col(x)

                    ## png(paste0('../plots/', prev.gene, '.heatmap.png'), height=8, width=8, units='in', res=150)
                    ## heatmap.2(x, trace='none', na.color='gray')
                    ## dev.off()

                    deg.fr <- ncol(x)
                    mean <- colMeans(x)
                    Sx <- cov(x)
                    D2 <- mahalanobis(x, mean, Sx, tol=1e-20)
                    ## pdf(paste0('../plots/', prev.gene, '.mahalanobis_d2.pdf'), height=8, width=8)
                    ## plot(density(D2, bw = 0.5), main=paste0('Mahalanobis distance, ', prev.gene))
                    ## rug(D2)
                    ## dev.off()

                    ## pdf(paste0('../plots/', prev.gene, '.qq.pdf'), height=8, width=8)		
                    ## qqplot(qchisq(ppoints(100), df = deg.fr), D2,
                    ##        main = expression("Qs-Q plot of Mahalanobis" * ~D^2 *
                    ##                          " vs. quantiles of" * ~ chi^2))
                    ## abline(0, 1, col = 'gray')
                    ## dev.off()

                    sig <- D2[D2>qchisq(0.999, ncol(x))]
                    if (length(sig) > 0) {
                            for (j in 1:length(sig)) {
                                    cat(names(sig)[j], prev.gene, log(pchisq(sig[j], df=deg.fr, lower.tail=F), 10), '\n', sep='\t')
                            }
                    }
                  }
		x <- NULL
	}
	prev.gene <- gene
	x <- rbind(x, v)
}





