#this is more of a note to myself:


#where x is a character array of genome positions, like chr1:1234-4321
do.parse<-function(x){unlist(strsplit(unlist(strsplit(x, ":")),"-"))}

posl <- lapply(x,do.parse)
bed <- do.call(rbind, posl)

write.table(bed, file="results.bed", sep="\t",quote=F)
