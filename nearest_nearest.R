#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(IRanges)

qw <- function(...) {
  as.character(sys.call()[-1])
}

args <- commandArgs(trailingOnly=TRUE)
resdir = args[1]

setwd(resdir)
res.tss<- read.csv("res_tss.csv")
res.mirna <- read.csv("res_mirna.csv")
res.exon <- read.csv("res_exon.csv")


#Create a filtered file of just the nearest for each peak.
ord<-order(res.tss[,"names"], abs(res.tss[,"tss.distancetoFeature"]), decreasing=FALSE)
res.tss.nearest<-res.tss[ord,]
res.tss.nearest<-res.tss.nearest[!duplicated(res.tss.nearest[,"names"]),]
write.csv(res.tss.nearest, file="res_tss_nearest.csv", row.names=F)


ord<-order(res.mirna[,"names"], abs(res.mirna[,"mirna.distancetoFeature"]), decreasing=FALSE)
res.mirna.nearest<-res.mirna[ord,]
res.mirna.nearest<-res.mirna.nearest[!duplicated(res.mirna.nearest[,"names"]),]
write.csv(res.mirna.nearest, file="res_mirna_nearest.csv", row.names=F)


ord<-order(res.exon[,"names"], abs(res.exon[,"exon.distancetoFeature"]), decreasing=FALSE)
res.exon.nearest<-res.exon[ord,]
res.exon.nearest<-res.exon.nearest[!duplicated(res.exon.nearest[,"names"]),]
write.csv(res.exon.nearest, file="res_exon_nearest.csv", row.names=F)
