library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("../Rcpp_funs.cpp") # Functions for computing the null distributions

# ===========================================================================

# Read in the simulation results (first run the script "scr.sh")
results <- read.table("data_props.txt", header=TRUE, fill=TRUE)
results <- results %>%
  mutate(RM_act = Muts - Seg_sites) %>%
  filter(!is.na(Pmin))

# Calculate the p-values for each row
N_iters <- 10000 # Number of iterations for computing each null distribution
Pen_factor <- 1.1 # Penalty factor F
P <- as.vector(t(read.table("mut_map.txt"))) # Mutation map (fitted to SARS-CoV-2 data)
m_range <- 250:450 # Range of number of sample segregating sites
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

# Calculate the p-values for the simulation results
pvals <- rep(NA, dim(results)[1])
for(i in 1:dim(results)[1]) {
  # Find the correct column number of "dists"
  ss <- which(dists[,52]==results$Seg_sites[i])
  rmfit <- results$Pmin[i]
  if(rmfit > 50) {
    rmfit <- 50
  }
  pvals[i] <- dists[ss, rmfit+1]
}
results <- cbind(results, pvals)

# Plot results (second panel of Figure 2)
results_plot <- results %>%
  group_by(Rec_rate) %>%
  summarise(reject = sum(pvals < 0.05)/n(),
            N = n()) %>%
  filter(Rec_rate > 0)

ggplot(data=results_plot, aes(x=Rec_rate, y=reject)) + 
  geom_point() +
  geom_line() +
  scale_x_continuous(name="Recombination rate", 
                     limits=c(0.0000001, 0.0000128), trans="log10") +
  scale_y_continuous(name="Proportion rejected", 
                     breaks=seq(0,1,0.1), labels=seq(0,1,0.1), limits=c(0, 1)) +
  theme_minimal() +
  annotation_logticks(sides="b")