#!/usr/bin/env Rscript

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

#these have already been mapped to their nearest transcript
peak.file <- args[1]
xpn.file <- args[2]
merge.file <- args[3]

peaks <- read.csv(peak.file)
xpn <- read.csv(xpn.file)


# the peaks are mapped to nearest ensembl gene IDs
# the expression data has Ensembl Gene ID annotation in the
# "Proportion_Ensembl_transcripts" field
ensembl.ids <- sub(".*(ENS.*)\\).*", "\\1", xpn[,"Proportion_Ensembl_transcripts"])
xpn<-cbind(xpn, ensembl.id=ensembl.ids)

# Chuck out anything that doesn't map to an Ensembl ID
peaks<-peaks[peaks[,"tss.feature"]!="",]
xpn<-xpn[xpn[,"ensembl.id"]!="",]

# Map between the 2 on the basis of Ensembl Gene ID.
merged <- merge(peaks, xpn, by.x="tss.feature", by.y="ensembl.id")

write.csv(merged, merge.file)




