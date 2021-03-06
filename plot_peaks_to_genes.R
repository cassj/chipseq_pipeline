#!/usr/bin/env Rscript

# This script takes a file as generated by, for example mm9RDtoGenes, containing
# peak data mapped to the nearest TSS and generates some summmary plots and information
# Output goes to the directory that contains the input file

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)

filename <-  args[1]

#for testing
filename <- "/mnt/esc/macs_300_1.0e-05/res_tss.csv"

peak.data <- read.csv(filename)

# > colnames(peak.data)
#  [1] "space"                        "start"                       
#  [3] "end"                          "width"                       
#  [5] "names"                        "values.Length"               
#  [7] "values.Summit"                "values.nTags"                
#  [9] "values.neg10log10pVal"        "values.FoldEnrichment"       
# [11] "values.FDR"                   "tss.strand"                  
# [13] "tss.feature"                  "tss.start_position"          
# [15] "tss.end_position"             "tss.insideFeature"           
# [17] "tss.distancetoFeature"        "tss.shortestDistance"        
# [19] "tss.fromOverlappingOrNearest" "mgi_symbol"                  
# [21] "description"  
