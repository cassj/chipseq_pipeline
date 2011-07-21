options(stringsAsFactors = FALSE);

args <- commandArgs(trailingOnly=TRUE)

#for testing
args <- c("results/peak_compare.csv", "TC", "/space/angela/Mash1_TC/results/alignment/bowtie/peakfinding/macs/macs_300_1.0e-05/res_tss_nearest.csv", "SC", "/space/angela/Mash1_SC/results/alignment/bowtie/peakfinding/macs/macs_300_1.0e-05/res_tss_nearest.csv")

comp <- read.csv(args[1])
a.name <- paste(args[2], ".csv", sep="")
a <- read.csv(args[3])
b.name <- paste(args[4], ".csv", sep="")
b <- read.csv(args[5])

a.merge <- merge(a, comp, by.x="names", by.y="id")
b.merge <- merge(b, comp, by.x="names", by.y="id") 

outdir <- dirname(args[1])

write.csv(a.merge, file=paste(outdir, a.name, sep="/"))
write.csv(b.merge, file=paste(outdir, b.name, sep="/"))
