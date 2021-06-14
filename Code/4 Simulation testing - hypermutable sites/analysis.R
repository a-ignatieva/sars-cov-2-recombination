# Code for performing simulation study to ascertain effect of adding hypermutable
# sites (Section 4.3.2, Figure S2, first panel of Figure 2)

library(tidyverse)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(latex2exp)
library(data.table)
library(smooth)
library(zoo)
library(R.utils)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../Rcpp_funs.cpp")

#=================================================================================

# The adjusted genomes are given in the files "mut_map_adjusted_p_Xf_Y.txt"
# where X = number of adjusted sites, and 
# Y = factor by which the mutation rate is multiplied.
# First, plot an example of the genome with added highly homoplasic sites (Figure S2)

# Mask problematic sites (see "1 Data processing, sampling and KwARG" for details of the masking)
vcf <- fread("../problematic_sites_sarsCov2.vcf", check.names=TRUE, skip="#CHROM")
homoplasic <- c(11083,13402,21575,16887,27384,3037,25563,8782,10323,11074,14408,6255,12880,21137,1059,1820,1912,6040,15324,29353,14805,15720,18060,28077,28826,28887,29253,29553,29700,6310,6312,11704,14786,17747,18756,20148,22661,23010,23403,29095,29422,3177,4084,6990,8078,11916,14724,14925,17247,18788,18877,20755,21648,24034,25947,26152,26461,27005,27046,27964,28881,29742,1457,4255,5784,7011,8293,8917,9223,10319,10507,11320,12781,13947,15760,16260,19684,22988,23422,24390,25916,26144,26530,26730,27525,28144,28311,28344,28851,28854,29751,379,541,833,884,1076,1570,1594,2113,3096,3253,3787,4113,4320,6573,7438,7765,10789,10851,11417,14747,15960,16762,17410,17639,17799,17858,18656,20031,20268,20275,21204,21707,23533,23587,24368,24389,24694,25494,25688,26211,26729,26735,28657,28688,28739,28857,28878,29540,29585,29734,313,490,1515,2455,2558,4809,6723,7479,8767,9477,9479,10097,10265,10450,11195,11801,13730,13929,14741,14912,15277,15927,16289,16381,17104,17373,17690,17944,18652,18713,18928,18998,19170,20931,23086,23707,23731,23929,24054,24862,25433,25572,25979,26124,26625,26936,27299,27635,27679,28580,28821,28836,28882,28883,29144,29635,29686)
mask <- c(vcf$POS, homoplasic)
mask <- unique(mask)

# Read in one of the mutation map files
P <- read.table("mut_map_adjusted_p_50f_2.txt")
position <- 1:29903
position <-  position[-mask]
mut_probs <- cbind(P, position)

# Figure S2 (without the colouring by base)
ggplot(data=mut_probs, aes(x=position, y=V1)) +
  geom_point(size=1) +
  theme_minimal() +
  xlab("Position") + 
  ylab(TeX("$\\widetilde{P}_{50,2}$")) +
  scale_x_continuous(breaks = seq(0,30000,5000), labels = seq(0,30000,5000), minor_breaks=seq(0,30000,1000)) +
  scale_y_continuous(breaks = seq(-0.00001,0.00014, 0.00002), labels = seq(-0.00001,0.00014, 0.00002), minor_breaks=NULL)

# =======================================================

# !! The file "data_props.txt" is created by running the script "scr.sh" (do this step now)

# Check the false positive rate vs number of added highly homoplasic sites
results <- read.table("data_props.txt", header=TRUE, fill=TRUE) %>%
  filter(!is.na(Pmin)) %>%
  mutate(RM_act = Muts - Seg_sites) %>%
  filter(RM_act >= 0) %>%
  mutate(Pmin = ifelse(Pmin > RM_act, RM_act, Pmin))

# Construct matrix giving the density of number of RMs from 0 to 50 (for a range of 
# number of sample segregating sites)
N_iters <- 10000 # Number of iterations for calculating each null distribution
Pen_factor <- 1.1 # penalty factor F
P <- read.table("mut_map_adjusted_p_0f_0.txt") # Unadjusted mutation map
m_range <- 250:450 # range of number of segregating sites to use
rm_max <- 51 

dists <- array(0, dim=c(length(m_range), rm_max))
for(i in 1:length(m_range)) {
  m <- m_range[i]
  s <- sim_null_dist(round(m*Pen_factor, 0), N_iters, P)
  h <- hist(s, breaks=seq(-0.5, max(s) + 0.5, 1))
  s <- h$density
  s <- c(s, rep(0, max(0, rm_max - length(s))))
  s <- rev(cumsum(rev(s)))
  dists[i,] <- s
}
dists <- cbind(dists, m_range)

# Calculating the p-values
pvals <- rep(NA, dim(results)[1])
for(i in 1:dim(results)[1]) {
  # Find the matching column number of dists
  ss <- which(dists[,52]==results$Seg_sites[i])
  rmfit <- results$Pmin[i]
  if(rmfit >= 50) {
    rmfit <- 50
  }
  pvals[i] <- dists[ss, rmfit+1]
}
results <- cbind(results, pvals)

plot_table <- results %>%
  group_by(j, f) %>%
  summarise(n_reject = sum(pvals < 0.05),
            N = n(),
            p_reject = n_reject/N)

# Results plot (first panel of Figure 2)
ggplot(data=plot_table, aes(x=as.factor(j), y=p_reject, 
                                  colour=as.factor(f), group=as.factor(f))) +
  geom_point(position=position_dodge(width=0.8), size=1.5) +
  geom_line(position=position_dodge(width=0.8), alpha=0.5) + 
  scale_y_continuous(breaks=seq(0,1,0.1), 
                     labels=seq(0,1,0.1), limits=c(0,1)) +
  theme_minimal() +
  xlab("Highly homoplasic sites") +
  ylab("Proportion rejected") +
  labs(colour = "H") +
  theme(legend.position = c(0.18,0.55), 
        legend.background = element_rect(fill="white", colour="white"))







