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

# Get list of sites to mask:
# https://github.com/W-L/ProblematicSites_SARS-CoV2
# https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473
# using update to 4 March 2021
vcf <- fread("problematic_sites_sarsCov2.vcf", check.names=TRUE, skip="#CHROM")
mask <- unique(vcf$POS)
# Also masking sites identified as homoplasic by van Dorp et al. (2020):
# https://doi.org/10.1016/j.meegid.2020.104351
# Supplementary table S5 (positions given in column 1)
homoplasic <- c(11083,13402,21575,16887,27384,3037,25563,8782,10323,11074,14408,6255,12880,21137,1059,1820,1912,6040,15324,29353,14805,15720,18060,28077,28826,28887,29253,29553,29700,6310,6312,11704,14786,17747,18756,20148,22661,23010,23403,29095,29422,3177,4084,6990,8078,11916,14724,14925,17247,18788,18877,20755,21648,24034,25947,26152,26461,27005,27046,27964,28881,29742,1457,4255,5784,7011,8293,8917,9223,10319,10507,11320,12781,13947,15760,16260,19684,22988,23422,24390,25916,26144,26530,26730,27525,28144,28311,28344,28851,28854,29751,379,541,833,884,1076,1570,1594,2113,3096,3253,3787,4113,4320,6573,7438,7765,10789,10851,11417,14747,15960,16762,17410,17639,17799,17858,18656,20031,20268,20275,21204,21707,23533,23587,24368,24389,24694,25494,25688,26211,26729,26735,28657,28688,28739,28857,28878,29540,29585,29734,313,490,1515,2455,2558,4809,6723,7479,8767,9477,9479,10097,10265,10450,11195,11801,13730,13929,14741,14912,15277,15927,16289,16381,17104,17373,17690,17944,18652,18713,18928,18998,19170,20931,23086,23707,23731,23929,24054,24862,25433,25572,25979,26124,26625,26936,27299,27635,27679,28580,28821,28836,28882,28883,29144,29635,29686)
# Plus additional sites to mask for SA data
mask <- c(mask, homoplasic, 22266:22745, 28254)
mask <- unique(mask)
mask <- sort(mask)

#===============================================================================

# Read in the the alignment
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


