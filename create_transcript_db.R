#!/usr/bin/env Rscript

# This script will create a TranscriptDB object for the specified species
# from the current version of Ensembl via Biomart and save it as an SQLite
# database in the specified output file

library("GenomicFeatures")

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

species <- args[1]
db.file <- args[2]

ensembl <- makeTranscriptDbFromBiomart(biomart="ensembl",
                                       dataset=paste(species, "_gene_ensembl", sep=""))

saveFeatures(ensembl, file=db.file)
