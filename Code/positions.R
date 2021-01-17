#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))

# Read in the VCF file containing SNP information
vcf <- fread(args[1], check.names=TRUE) %>%
  select(POS)
write.table(t(vcf), file=args[2], sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)











