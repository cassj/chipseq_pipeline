###
#
# Script to compare 2 ChIPseq bam files to identify significantly different peaks
#
###

# you will need:
# The BAM files for each of the samples
# A BED file containing the peak regions for each of the samples


#!/usr/local/bin/Rscript compare_peaks /path/to/outdir N_threads /path/to/bam1 /path/to/bed1 /path/to/bam2 /path/to/bed2 ...

options(stringsAsFactors = FALSE);


args <- commandArgs(trailingOnly=TRUE)

#for testing
args<-c("/mnt/data/",4, "/mnt/astro/BAM/IP.bam" ,"/mnt/astro/Macs/NA_peaks.bed", "/mnt/esc/chip_export_sorted_nodups.bam","/mnt/esc/macs_300_1.0e-05/EscChIPseqREST_peaks.bed")

resdir <- args[1]
threads <- args[2]
args <- args[-1:-2]

inds <- 1:length(args)
bam.files <- args[which(inds%%2!=0)]
bed.files <- args[which(inds%%2==0)]


library(IRanges)
library(ShortRead)
library(snow)
library(baySeq)
library(DESeq)
library(rtracklayer)
library(Rsamtools)

# For ChIPseq data we aren't dealing with that many locations, so probably
# we're good using just the number of cores on the AWS machine. 
if(is.null(threads) || threads==1){
  cl <- NULL
}else{
  cl <- makeCluster(threads,"SOCK")
}


# Read in theBED files as RangedData.
# Retrieve the read data for those regions
# Build the counts table
# This seems to take about an hour for 4K regions.

beds <- list()
counts <- NULL
seglens <- NULL
for(i in 1:length(bed.files)){
  beds[[i]] <- import(bed.files[i])

  #get the count data for these ranges from each bam file
  bam.counts <- list()
  for(j in 1:length(bam.files)){
    what <- c("qname") 
    param <- ScanBamParam(which=beds[[i]], what=what)
    bam <- scanBam(bam.files[[j]], param=param)
    bam.counts[[j]] <- sapply(bam, function(x){length(x$qname)})
  }

  beds[[i]] <- as.data.frame(beds[[i]])
  nms <- paste(beds[[i]][,"space"], paste(beds[[i]][,"start"], beds[[i]][,"end"], sep="-"), sep=":")
  these <- do.call(cbind, bam.counts)
  colnames(these) <- bam.files
  if(is.null(counts)){
    counts <- these
  }else{
    counts <- rbind(counts, these)
  }
  if(is.null(seglens)){
    seglens <- beds[[i]][,"end"]-beds[[i]][,"start"]+1
  }else{
    seglens <- c(seglens,beds[[i]][,"end"]-beds[[i]][,"start"]+1)
  }
}

save(counts, file=paste(resdir,"counts.RData", sep="/"))

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
NBML <- getPriors.NB(cd, samplesize=1000, estimation="QL", cl=cl)
post.NBML <- getLikelihoods.NB(NBML,pET="BIC",cl=cl)

props <- post.NBML@estProps

#very small fraction are different. This is < 1.
#Which either means that we haven't got enough data or there's no difference.
props[2]*nrow(counts)

#Try it on the Mash data maybe?
