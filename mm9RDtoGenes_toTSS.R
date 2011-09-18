#!/usr/local/bin/Rscript

options(stringsAsFactors = FALSE);

library(biomaRt)


library(IRanges)

library(ChIPpeakAnno)

qw <- function(...) {
  as.character(sys.call()[-1])
}


args <- commandArgs(trailingOnly=TRUE)
filename = args[1]

##for testing
#filename = "peaks.RangedData.RData"

rd <- get(load(filename))


ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

#remove "chr" prefix for  ChIPpeakAnno
names(rd)<-gsub("chr","",names(rd))

#change "M" to "MT" for ChIPpeakAnno
id<-which(names(rd)=="M")
if (length(id)>0){
   names(rd)[id]<-"MT"
}


# NOTE: TSS.mouse.NCBIM37 is actually *gene* start and end positions, not individual transcripts.

#get the most recent annotation data from ensembl
tss <- getAnnotation(ensmart, "TSS")
save(tss, file=paste(dirname(filename),"/tss.RData",sep=""))

#chippeakanno loses peak values - merge them back on later
rd.df <- as.data.frame(rd)
vals <- rd.df[,c("names","values.nTags","values.neg10log10pVal","values.FoldEnrichment","values.FDR")]
colnames(vals) <- c("Names","nTags","neg10log10pVal","FoldEnrichment","FDR")

#find peaks that are near to genes, measure distance to TSS (feature start on +ve and end on -ve strands), return all near peaks and overlapping - hope that it does the strand right....this will need testing.

peak_to_gene <- annotatePeakInBatch(rd,
                                    AnnotationData=tss,
                                   PeakLocForDistance = "middle",    # from the middle of the peak
                                    FeatureLocForDistance = "TSS",  # to the TSS of the gene
                                    output = "both",
                                    multiple=TRUE
                                    )

#save as takes ages to make
save(peak_to_gene, file = paste(dirname(filename),"/peak_to_gene.RData",sep=""))
#peak_to_gene <- get(load("peak_to_gene.RData"))

#convert to dataframes

peak_to_gene.df <- as.data.frame(peak_to_gene)

#annotate genes from ensembl IDs using BiomaRt

filters <- c("ensembl_gene_id")
values<-unique(peak_to_gene.df[,"feature"])
attributes <- c("ensembl_gene_id","mgi_symbol", "description")

annot <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)


#ditch any that don't have a symbol or a description
no.anno <- intersect(
                     which(annot[,"mgi_symbol"]==""),
                     which(annot[,"description"]==""))
annot <- annot[-1*no.anno,]


# a few have multiple bits of annotation
annot<-cbind(annot, alt.annot="")
dups<-annot[duplicated(annot[,"ensembl_gene_id"]), "ensembl_gene_id"]

#keep the first one and add all the others as alt.annot 
for (d in dups){
  inds <- which(annot[,"ensembl_gene_id"]==d)
  this.alt.annot <- annot[inds[-1], c("mgi_symbol", "description")]
  annot[inds[1],"alt.annot"] <- paste(paste(this.alt.annot[,1], this.alt.annot[,2]), collapse="; ")
}
annot <- annot[!duplicated( annot[,"ensembl_gene_id"] ), ]
rownames(annot) <- annot[,"ensembl_gene_id"]

#merge peaks back to genes by Ensembl IDs

merged.res <- merge(peak_to_gene.df, annot, by.x = "feature", by.y = "ensembl_gene_id")

#merge back to peak values by peak name

merged2.res <- merge(merged.res,vals, by.x = "peak", by.y = "Names")

#rearrange, save and test
res <- merged2.res[,qw(peak,space,start,end,width,nTags,neg10log10pVal,FoldEnrichment,FDR,feature,strand,start_position,end_position,insideFeature,distancetoFeature,shortestDistance,fromOverlappingOrNearest,mgi_symbol,description)]
colnames(res) <- qw(Peak,Chromosome,Peak_start,Peak_end,Peak_width,nTags,neg10log10pVal,FoldEnrichment,FDR,EnsemblID,Gene_strand,Gene_start,Gene_end,insideFeature,distancetoFeature,shortestDistance,fromOverlappingOrNearest,Symbol,Description)

res <- res[order(res[,"neg10log10pVal"],decreasing=TRUE),]
write.csv(res, file = paste(dirname(filename),"/nearest_or_overlapping_peak_to_gene.csv",sep=""))

#find nearest peak to each gene - take nearest peak only

res.o <- res[order(abs(res[,"distancetoFeature"]),decreasing=FALSE),]

res.od <- res.o[which(!duplicated(res.o[,"Peak"])),]

res.odd <- res.od[which(!duplicated(res.od[,"EnsemblID"])),]

write.csv(res.odd, file = paste(dirname(filename),"/nearest_peak_to_gene_TSS.csv",sep=""))













