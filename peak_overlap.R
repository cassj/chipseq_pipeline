# I still maintain that trying to do this as a naive overlap
# is misleading, given that you're forcing a dichotomous structure on
# continuous data, but meh.

# We also need to be sure that the cut-off lists have been generated using the same FDR. Should probably do some kind of sanity check here I guess.

library(IRanges)

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

output.dir = args[1]
args <- args[-1]



#for testing
output.dir <- "results"
args <- c(
          "Astro", "/space/cassj/REST_ChIP_mouse_astrocytes_volta/Macs/results/NA_peaks.RangedData.RData",
          "ESC", "/space/cassj/mouse_ESC_REST_ChIPseq/results/alignment/bowtie/peakfinding/macs/macs_300_1.0e-05/EscChIPseqREST_peaks.RangedData.RData",
          "NSC", "/space/cassj/mouse_NSC_REST_ChIP/results/alignment/bowtie/peakfinding/macs/macs_300_1.0e-05/NscChIPseqREST_peaks.RangedData.RData"
          )

filenames <- args[which( (1:length(args) %% 2) == 0)]
nms <- args[which( (1:length(args) %% 2) == 1)]


rds <- lapply(filenames, function(x){get(load(x))})
names(rds) <- nms



# these are the combinations we need to calculate overlaps for
combos <- lapply(2:length(nms), function(x){combn(nms,x)})

# these are the 2-way combinations, which we calculate with findOverlaps

# overlap 2 ranged data objects and convert the resulting
# RangesMatchingList into a dataframe we can run merge on.
# Add the actual ranges to the dataframe too so we can map
# back to annotations etc.
do.overlaps <- function(x,y, nm.x, nm.y){
   rml <- findOverlaps(x,y)
   non.empty <- sapply(rml, function(x){nrow(x@matchMatrix)})!=0
   rml <- rml[non.empty]
   mms <- lapply(names(rml), function(m){
     mm <-  rml[[m]]@matchMatrix
     data.frame(chr=m,
                mm,
                as.data.frame(ranges(x[m][mm[,1],]))[,"names"],
                as.data.frame(ranges(y[m][mm[,2],]))[,"names"]
                )  })
   res <- do.call(rbind, mms)
   colnames(res) <- c("chr", paste(c(nm.x, nm.y),"index",sep="."), paste(c(nm.x, nm.y),"pos",sep="."))

   # get the actual regions from the original objects
   
   return(res)
   
}

# overlaps names the cols "query" and "subject, so we need to pass the RD objects and names for them
ols <- list(apply(combos[[1]], 2, function(x){do.overlaps(rds[[x[1]]],rds[[x[2]]], x[1], x[2]) } ))




# now use the 2-way overlaps to construct the higher level combinations
for(i in 2:length(combos)){
  these <- combos[[i]]
  these.ols <- list()
  for(j in 1:ncol(these)){
    #work out which 2 of the previous row you need to merge to get this one
    #eg, if ABCD, we need ABC and BCD from the previous row:
    this <- these[,j]
    first.half <- this[1:(length(these)-1)]
    second.half <- this[2:length(these)]
    x.ind <- which(apply(combos[[i-1]],2,function(x){all(x==first.half)}))
    y.ind <- which(apply(combos[[i-1]],2,function(x){all(x==second.half)}))
    x <- ols[[i-1]][[x.ind]]
    y <- ols[[i-1]][[y.ind]]
    these.ols[[j]] <- merge(x,y)
  }
  ols[[i]] <- these.ols

}


#save the overlap data just in case
save(ols, file=paste(output.dir, "ols.RData", sep="/"))


# Venn diagram table:

# We need to give each region of overlap a single name
# Start from the overlap of everything and work back

n <- length(ols)
pos.nms <- ols[[n]][[1]]
pos.nms <- pos.nms[,grep('.pos',colnames(pos.nms))]

make.name <- function(x){
  splt <- unlist(strsplit( unlist(strsplit(as.character(x),':')), '-'))
  inds <- grep("chr", splt)
  chr  <- splt[inds[1]]
  mn <- min(as.numeric(splt[-1*inds]))
  mx <- max(as.numeric(splt[-1*inds]))
  nm <- paste(chr, paste(mn,mx, sep="-"), sep=":")
}



# *** this isn't working. Gnorks. *** 
# generate a name lookup table
lookup <- list()
for (i in n:1){   # 1..n levels in ols (indexed by i)
  m <- length(ols[[i]]) 
  for(j in 1:m){  # 1..m overlaps at this level (indexed by j)
    these <- ols[[i]][[j]]
    these <- these[,-1*grep(".index", colnames(these))]
    these.nms <- apply( these, 1, make.name )
    for(r in 1:nrow(these)){  #for each row (r) in this overlap 
      for(p in as.character(these[r,grep(".pos", colnames(these))]) ){ # and for each position (p) in this row  
        if(is.null(lookup[[get("p")]])){
          lookup[[get("p")]]<- these.nms[r]
        }
      }
    }
  }
}

lookup<-unlist(lookup)
ol.names <- unique(lookup)


bool.table<-matrix(nrow=length(ol.names), row.names=ol.names, ncol=


#venn.df <- data.frame(


