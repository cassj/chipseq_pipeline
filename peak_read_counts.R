#!/usr/bin/env Rscript

###
#
# For each region in the BED files, fetches the read count in that regions from the BAM file.
#
###

# you will need:
# The BAM files for each of the samples
# A BED file containing the peak regions for each of the samples


#!/usr/local/bin/Rscript peak_read_counts.R /path/to/outdir N_threads /path/to/bam1 /path/to/bed1 /path/to/bam2 /path/to/bed2 ...

options(stringsAsFactors = FALSE);


args <- commandArgs(trailingOnly=TRUE)

#for testing
#args<-c("/mnt/data/",4, "/mnt/astro/BAM/IP.bam" ,"/mnt/astro/Macs/NA_peaks.bed", "/mnt/esc/chip_export_sorted_nodups.bam","/mnt/esc/macs_300_1.0e-05/EscChIPseqREST_peaks.bed")
#args <- c("/mnt/data/",4,"/mnt/TC/CMN066_s_8_export_sorted_nodups.bam", "/mnt/TC/macs_300_1.0e-05/Mash1_TC_peaks.bed", "/mnt/SC/chip_export_sorted_nodups.bam", "/mnt/SC/macs_300_1.0e-05/Mash1_SC_peaks.bed")


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


# Read in the BED files as RangedData.
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
  col.nms <- sub("-","", gsub("/","-", bam.files))
  colnames(these) <- col.nms
  these <- data.frame(sample=bed.files[i], these)
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



# We'll need the library sizes later, so get them from the summary files
# Note that these should be generated using the samtools flagstat program
summaries <- gsub('.bam', '.summary', bam.files )
libsizes <- sapply(summaries, function(x){
  l <- readLines(x,1)
  l <- sub("\\s+.*","",l, perl=T)
})


save(libsizes, file=paste(resdir,"libsizes.RData", sep="/"))
