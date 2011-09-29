#!/usr/bin/env Rscript

# This script takes a BAM file and a TranscriptDB sqlite file and plots
# the distribution of reads across an "average" gene.

library("GenomicFeatures")
library("GenomicRanges")
library("Rsamtools")

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

bam <- args[1]
transcripts <- args[2]

#for testing
bam <- "/mnt/esc/chip_export_sorted_nodups.bam"
db.file <- "/mnt/data/mmusculus_ensembl_20110927.sqlite"

genes <- loadFeatures(db.file)
genes <-as.list(genes)


# load the bam into a GappedAlignments object ?
aln <- readGappedAlignments(bam)

# aln has chr, exonList doesn't.
rname(aln) <- gsub("chr","",rname(aln))
# and uses MT not M
rname(aln) <- gsub("M","MT",rname(aln))


# This will give you an RLE along each chr. 
cov <- coverage(aln)

get.cov <- function(x){
  start <- as.integer(x["tx_start"]) - 20000
  end <- as.integer(x["tx_end"]) + 20000
  chr <- x["tx_chrom"]
  return(cov[[as.character(chr)]][start:end])
}

test<-apply(genes[[1]], 1, get.cov)



what <-  c("rname", "pos", "length")
which <- RangesList(space=genes[[1]][,"tx_chrom"], start=genes[[1]][,"tx_start"], end=genes[[1]][,"tx_end"])


  
# This doesn't do what they think it does.
counter <- function(gnModel)
{
    hits <- countOverlaps(aln, gnModel) 
    counts <- countOverlaps(gnModel, aln[hits==1])
    names(counts) <- names(gnModel)
    counts
}





counts <- counter(exonList[[1]])
save(counts, file=file.path(outputDir, "counts.rda"))


### faffing about #####

# get transcripts
transcripts <- transcripts(genes)

# get transcripts by gene
transcriptList <- transcriptsBy(genes, by="gene")

# get exons
exons <- exons(genes)

# get exons by transcripts
exonList <-  exonsBy(genes, by="gene")

# get introns by genei
ntronList <- intronsByTranscript(genes)


transcriptsByOverlap
exonsByOverlap




