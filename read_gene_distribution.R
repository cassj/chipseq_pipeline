#!/usr/bin/env Rscript

# This script takes a BAM file and a TranscriptDB object and plots the distribution.


options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

bam <- args[1]

#for testing
bam <- "/mnt/esc/chip_export_sorted_nodups.bam"
