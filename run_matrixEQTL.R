#!/usr/bin/env Rscript

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.


# print usage
usage <- function() {
    cat(
          'usage: run_matrixEQTL.R <cis distance> <tissue> <directory> <out>
run_matrixEQTL.R
author: Colby Chiang (colbychiang@wustl.edu)
    cis distance cis distance
    tissue       tissue
    directory    directory for input and output
    out          prefix for output file
')
}

# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
CIS_DISTANCE <- as.numeric(args[1])
TISSUE <- args[2]
DIR <- args[3]
OUT_PREFIX <- args[4]

# Check input args
if (is.na(CIS_DISTANCE)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(TISSUE)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(DIR)) {
  usage()
  quit(save='no', status=1)
}
if (is.na(OUT_PREFIX)) {
  usage()
  quit(save='no', status=1)
}
    

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
# install.packages("MatrixEQTL")
## install.packages('/gscmnt/gc2719/halllab/users/cchiang/src/MatrixEQTL', repos = NULL, type="source")
library('MatrixEQTL')

## Location of the package with the data files.
base.dir = DIR;
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0(base.dir, "/", TISSUE, ".sv.scaled.sup_10_samp.txt");
snps_location_file_name = paste0(base.dir, "/", TISSUE, ".varsloc.txt.gz");

# Gene expression file name
expression_file_name = paste0(base.dir, "/", TISSUE, ".expr.txt");
gene_location_file_name = "/gscmnt/gc2719/halllab/users/cchiang/projects/gtex/data/GTEx_Analysis_2015-01-12/eqtl_data/eQTLInputFiles_genePositions/GTEx_Analysis_2015-01-12_eQTLInputFiles_genePositions.txt.gz"
## gene_location_file_name = paste0(base.dir, "/", TISSUE, ".genesloc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir, "/", TISSUE, ".covariates.txt");

# Output file name
output_file_name_cis = paste0(OUT_PREFIX, ".cis_eqtl.txt");
output_file_name_tra = paste0(OUT_PREFIX, ".trans_eqtl.txt");

# Only associations significant at this level will be saved
## pvOutputThreshold_cis = 1;
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-5;
## pvOutputThreshold_tra = 1e-3;


# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
## cisDist = 1e6;
cisDist = CIS_DISTANCE;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


# detach('package:MatrixEQTL', unload=TRUE)
# install.packages('/gscmnt/gc2719/halllab/users/cchiang/src/MatrixEQTL', repos = NULL, type="source")
# library('MatrixEQTL')

me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);

# min.pv.by.genesnp = TRUE

## unlink(output_file_name_tra);
## unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
## show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
## show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me)
