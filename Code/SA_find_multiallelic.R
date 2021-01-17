#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(data.table)))

#===============================================================================

# Get list of sites to mask from Virological:
# https://github.com/W-L/ProblematicSites_SARS-CoV2
# https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
vcf <- fread("problematic_sites_sarsCov2.vcf", check.names=TRUE, skip="#CHROM")
mask <- unique(vcf$POS)
mask <- c(mask, 22339:22523, 22287)
mask <- unique(mask)

#===============================================================================

# Read in the multiple alignment
dat <- read.table(args[1], header=FALSE, fill=TRUE, 
                  row.names=NULL, sep="", strip.white=TRUE,
                  colClasses = "character")

# Split names from data
L <- dim(dat)[1]
dat_names <- dat[seq(1,L,2),] 
dat_names$N = seq(1,L,2)
dat_values <- dat[seq(2,L,2),]
dat_values$N = seq(2,L,2)

# Get multi-allelic sites and count number of SNPs
a <- c()
ct_mult <- 0
ct_snps <- 0
ct_snps_masked <- 0
mult_sites <- c()
snp_sites <- c()
for(i in 1:(dim(dat_values)[2]-1)) {
  lev = dat_values[,i] %>% unique()
  lev <- lev[lev %in% c("A", "C", "T", "G")]
  c = length(lev)
  if(i %in% mask) {
    dat_values[,i] <- "N"
    if(c == 2) {
      ct_snps_masked <- ct_snps_masked + 1
    }
  }
  else {
  if(c == 2) {
    ct_snps <- ct_snps + 1
    snp_sites <- c(snp_sites, i)
  }
  if(c > 2) {
    a <- c(a, i)
    ct_mult <- ct_mult + 1
    mult_sites <- c(mult_sites, i)
  }
  }
}
print(paste0("SNP sites: ", ct_snps), quote = FALSE)
print(paste0("Plus masked: ", ct_snps_masked), quote = FALSE)
print(paste0("Plus multi-allelic: ", ct_mult), quote = FALSE)
print("SNP sites:", quote = FALSE)
snp_sites
print("Multi-allelic sites:", quote = FALSE)
mult_sites
print("Masked sites:", quote = FALSE)
mask
write.table(t(snp_sites), file=args[3], sep=" ", na="", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Mask the multi-allelic sites
for(i in mult_sites) {
  dat_values[,i] <- "N"
}

# Combine with the sequence names
dat2 <- rbind(dat_names, dat_values)
dat2 <- dat2[order(dat2$N),]
dat2 <- subset(dat2, select = -c(N))

# Write to file
write.table(dat2, file=args[2], sep="", na="", quote=FALSE, row.names=FALSE, col.names=FALSE)


