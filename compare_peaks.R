###
#
# Script to compare 2 ChIPseq bam files to identify significantly different peaks
#
###

# you will need:
# The BAM files for each of the samples
# A BED file containing the peak regions for each of the samples


#!/usr/local/bin/Rscript compare_peaks /path/to/counts_file.RData /path/to/libsizes.RData

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

library(IRanges)
library(ShortRead)
library(snow)
library(baySeq)
library(DESeq)
library(rtracklayer)
library(Rsamtools)
library(edgeR)

counts.file <- args[1]
libsizes.file <- args[2]

#for testing
counts.file <- "/mnt/data/counts.RData"
libsizes.file <- "/mnt/data/libsizes.RData"

counts <- get(load(counts.file))
libsizes <- get(load(libsizes.file))
libsizes <- as.integer(libsizes)

# assuming no reps at the moment
groups <- 1:(ncol(counts) - 1)


### edgeR

# I think, as we have no reps, we get common dispersion of 0 which is essentialy Poisson.
# Which means the p values are far too low cos the Poisson doesn't model it very well.
d <- DGEList(counts[,-1], group=groups)

# Can we use the input sample to generate a null dist?

#unsure if we need to worry about normalisation. Supposed to correct for highly expressed transcripts. Don't think we're dealing with the same range in chip
#Ignore this step for now


#Estimating dispersion for the NB distribution
#quantile adjusted conditional maximum likelihood (qCML) for single factors
#Cox-Reid Profile-adjusted likelihood (CR) for multiple factors


# DE
# exact test -> based on qCML
# GLM

### DESeq


counts.table <- counts[,-1]
group.names <- c("TC","SC")

counts.table <- counts[,-1]
counts.table <- counts.table + 1
cds <- newCountDataSet( counts.table, group.names)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds, method="blind")
res <- nbinomTest(cds, "TC","SC")

res <- res[order(res[,"pval"]),]


# need to map res back to genes, but it kinda works.
# Really not getting enough sig with no reps though.


# ok, this kinda works but we're getting infinite fold changes. How to deal with the 0 values?

### baySeq

# can't get this to work with multiple threads
if(is.null(threads) || threads==1){
  cl <- NULL
}else{
  cl <- makeCluster(threads,"SOCK")
}



# Now use baySeq to determine differentially expressed counts.
replicates <- c(1,2)
groups <- list(NDE=c(1,1), DE=c(1,2))

# we need to get the library sizes from the summary files
summaries <- gsub('.bam', '.summary', bam.files )
libsizes <- sapply(summaries, function(x){
  l <- readLines(x,1)
  l <- sub("\\s+.*","",l, perl=T)
})



cd <- new("countData", data=counts, replicates=replicates, libsizes=as.integer(libsizes), groups=groups, seglens=seglens)

# Determine significance of differential counts assuming Negative Binomial count distribution.

# Bootstrap an empirical distribution from the data.
# This should really use 10,000 iterations, Maybe best to use the multiple cores
NBML <- getPriors.NB(cd, samplesize=10000, estimation="QL", cl=cl)
post.NBML <- getLikelihoods.NB(NBML,pET="BIC",cl=cl)

props <- post.NBML@estProps

#very small fraction are different. This is < 1.
#Which either means that we haven't got enough data or there's no difference.
props[2]*nrow(counts)

#Try it on the Mash data maybe?
