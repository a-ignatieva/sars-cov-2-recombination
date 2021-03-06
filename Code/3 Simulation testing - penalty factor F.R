# Code for performing simulation study in Section S4.3.1 and producing Figure S1

library(tidyverse)
library(data.table)
library(zoo)
library(wavethresh)
library(EbayesThresh)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("Rcpp_funs.cpp") # Functions for simulating the null distribution

# ===============================================================

# Function to simulate P as an autoregressive process
sim_P <- function(L=29903, sigma=0.000001, m=100, A, C, G, U) {
  P <- arima.sim(model=list(order=c(1,0,1), ar = 0.99999, ma = 0.001, sd = sigma), n = L)
  P <- P - min(P)
  P <- P/sum(P)
  
  # Correct the proportions of sites undergoing mutation to match SARS-CoV-2
  # (to add another source of rate heterogeneity)
  P[A] <- 0.25 * P[A]/sum(P[A])
  P[C] <- 0.27 * P[C]/sum(P[C])
  P[G] <- 0.23 * P[G]/sum(P[G])
  P[U] <- 0.25 * P[U]/sum(P[U])
  
  return(P)
}

# ===============================================================

# Get a list of positions by base type in the reference sequences
snp <- data.frame(base = t(read.table("ref_test.txt", header=FALSE, fill=TRUE, row.names=NULL, sep = " ", colClasses="character"))) %>%
  mutate(position = 1:29903)

positions_list <- list("A" = snp$position[snp$base == "A"], 
                       "C" = snp$position[snp$base == "C"], 
                       "G" = snp$position[snp$base == "G"], 
                       "T" = snp$position[snp$base == "T"])

# ===============================================================

mutations <- 20000 # number of mutations for generating the dataset
n <- 1000 # number of reps for null dist simulation
N <- 10 # number of reps for each run
m_range <- c(100, 300, 500) # number of sample segregating sites
p_range <- c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5) # penalty factor F
plotting_freq = 0 # how often the results should be plotted while the simulation is running

stats <- array(NA, dim=c(12, N*length(m_range)*length(p_range)))
# Will store the results in a table. The rows are:
# Number of mutations, number of seg sites, 95% quantile (true), 95% quantile (fitted)
# 90% quantile (true), 90% quantile (fitted), mean (true), mean(fitted), median (true), 
# median (fitted), sample number of seg sites
counter <- -1
for(k in 1:length(p_range)) {
  p <- p_range[k]
  for(j in 1:length(m_range)) {
    m <- m_range[j]
    counter <- counter + 1
    for(i in (1:N)+counter*N) {
      
      if(plotting_freq != 0 & i%%plotting_freq == 0) {
        plotting <- TRUE
        pdf(file=paste0("plot_", i, "_m_", m, ".pdf"), width=10, height=5)
      }
      else {
        plotting <- FALSE
      }
      
      # Simulate a vector P, with mutation rate varying by base and along the genome:
      Ptrue <- sim_P(L=29903, sigma = 0.000001, m = 1000, 
                     A=positions_list[["A"]], 
                     C=positions_list[["C"]], 
                     G=positions_list[["G"]], 
                     U=positions_list[["T"]])
      
      # Simulate some data using this true vector of probabilities P
      sites <- sample(1:29903, size=mutations, replace=TRUE, prob=Ptrue)
      stats[1, i] <- length(sites)
      tab <- table(sites)
      sites <- unique(sites)
      stats[2, i] <- length(sites)
      data <- rep(0, 29903)
      data[sites] <- 1
      
      # Perform the wavelet decomposition to fit estimate of mutation map
      mut_probs <- data.frame(position = 1:29903, 
                              rate = rep(NA, 29903), 
                              base = as.character(snp$base), 
                              stringsAsFactors=FALSE)
      for(b in c("A", "C", "G", "T")) {
        X <- data[positions_list[[b]]]
        Xpos <- positions_list[[b]]
        Xpos <- Xpos[!is.na(X)]
        X <- X[!is.na(X)]
        
        # Padding to make total length a power of 2
        # This is done by reflecting the data at the endpoints
        pow <- ceiling(log(length(X), base=2))
        ext_left <- floor((2^pow - length(X))/2)
        ext_right <- length(X) - ceiling((2^pow - length(X))/2) + 1
        Xf <- c(rev(X[1:ext_left]), X, rev(X[length(X):ext_right]))
        
        # Wavelets with empirical Bayes thresholding
        Y <- wd(Xf, filter.number=6, family="DaubLeAsymm")
        Y <- ebayesthresh.wavelet(Y, threshrule="soft")
        Y <- wr(Y)
        Y <- Y[(ext_left+1):(ext_left + length(X))]
        
        mut_probs$rate[mut_probs$position %in% Xpos] <- Y
      }
      
      mut_probs <- mut_probs %>%
        mutate(rate = ifelse(rate < 0, 0, rate),
               rate_norm = rate/sum(rate, na.rm=TRUE))
      Y <- mut_probs$rate_norm
      
      # Plot true vs fitted P vector
      if(plotting) {
        par(mfrow=c(2,3))
        plot(Ptrue, type='l', col=alpha("blue", 0.1), ylab="Mut prob", xlab="Position")
        points(mut_probs$rate_norm, type='l', col=alpha("red", 0.1))
   
        plot(Ptrue[10000:10300], type='p', 
             ylim=c(0, max(c(Ptrue[10000:10300], Y[10000:10300]))), 
             col=snp$base[10000:10300])
        plot(Y[10000:10300], type='p', 
             ylim=c(0, max(c(Ptrue[10000:10300], Y[10000:10300]))), 
             col=snp$base[10000:10300])
        
        plot(sort(Ptrue), type='l')
        points(sort(Y), type='l', col="red")
      }
      
      # Simulate null distribution using Ptrue and the fitted estimate
      p1 <- sim_null_dist(m, n, Ptrue)
      p2 <- sim_null_dist(ceiling(m*p), n, Y)
      
      if(plotting) {
        h1 <- hist(p1, freq=FALSE, breaks=seq(-0.5, max(c(p1,p2))+0.5, 1))
        h2 <- hist(p2, freq=FALSE, breaks=seq(-0.5, max(c(p1, p2))+0.5, 1))
        dev.off()
      }
      
      # Save the results
      stats[3, i] <- quantile(p1, 0.95)
      stats[4, i] <- quantile(p2, 0.95)
      stats[5, i] <- quantile(p1, 0.90)
      stats[6, i] <- quantile(p2, 0.90)
      stats[7, i] <- mean(p1)
      stats[8, i] <- mean(p2)
      stats[9, i] <- median(p1)
      stats[10, i] <- median(p2)
      stats[11, i] <- m
      stats[12, i] <- p
      
      print(i)
    }
  }
}

# ===============================================================
# Figure S1
# ===============================================================

validity <- as.data.frame(t(stats))
colnames(validity) <- c("total_muts", "seg_sites", "q95_true", "q95_fit", "q90_true", "q90_fit", "mean_true", "mean_fit", "median_true", "median_fit", "sample_sites", "slack_prop")

# Calculate differences between true and fitted distributions
validity <- validity %>%
  filter(!is.na(total_muts)) %>%
  group_by(sample_sites, slack_prop) %>%
  mutate(q95_true = ceiling(q95_true),
         q95_fit = ceiling(q95_fit),
         q90_true = ceiling(q90_true),
         q90_fit = ceiling(q90_fit),
         median_true = ceiling(median_true),
         median_fit = ceiling(median_fit),
         mean_diff = mean_true - mean_fit,
         median_diff = median_true - median_fit,
         q95_diff = q95_true - q95_fit,
         q90_diff = q90_true - q90_fit,
         q95_fails = as.numeric(q95_true > q95_fit),
         q90_fails = as.numeric(q90_true > q90_fit),
         median_fails = as.numeric(median_true > median_fit)) %>%
  ungroup()

counts <- validity %>%
  group_by(sample_sites, slack_prop) %>%
  summarise(q90_fails = sum(q90_fails)/500,
            q95_fails = sum(q95_fails)/500,
            median_fails = sum(median_fails)/500,
            q90_fails = ifelse(q90_fails > 0, q90_fails, NA),
            q95_fails = ifelse(q95_fails > 0, q95_fails, NA),
            median_fails = ifelse(median_fails > 0, median_fails, NA),
            q90_y_pos = min(q90_diff)-1.5,
            q95_y_pos = min(q95_diff)-1.5,
            median_y_pos = min(median_diff)-1.5)

# Figure S1 (middle)
ggplot(data=validity, aes(x=as.factor(slack_prop), y=q90_diff, group=as.factor(sample_sites))) +
  geom_count(aes(colour=as.factor(sample_sites)), position=position_dodge(width=1), show.legend=FALSE) +
  xlab("Penalty factor") +
  ylab("Difference with true 90th percentile value") +
  labs(colour = "Sample\nseg sites\n", size = "Count") +
  theme_minimal() +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_text(data=counts, aes(x=as.factor(slack_prop), 
                             label=scales::percent(q90_fails, accuracy=0.1), 
                             colour=as.factor(sample_sites), ymax=3, y = q90_y_pos), 
            angle=-90, 
            position = position_dodge(width = 1), show.legend = FALSE)

# Figure S1 (right)
ggplot(data=validity, aes(x=as.factor(slack_prop), y=q95_diff, group=as.factor(sample_sites))) +
  geom_count(aes(colour=as.factor(sample_sites)), position=position_dodge(width=1)) +
  xlab("Penalty factor") +
  ylab("Difference with true 95th percentile value") +
  labs(colour = "Sample\nseg sites\n", size = "Count") +
  theme_minimal() +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_text(data=counts, aes(x=as.factor(slack_prop), 
                             label=scales::percent(q95_fails, accuracy=0.1), 
                             colour=as.factor(sample_sites), ymax=3, y = q90_y_pos), 
            angle=-90, 
            position = position_dodge(width = 1), show.legend = FALSE)

# Figure S1 (left)
ggplot(data=validity, aes(x=as.factor(slack_prop), y=median_diff, group=as.factor(sample_sites))) +
  geom_count(aes(colour=as.factor(sample_sites)), position=position_dodge(width=1), show.legend=FALSE) +
  xlab("Penalty factor") +
  ylab("Difference with true median") +
  labs(colour = "Sample\nseg sites\n", size = "Count") +
  theme_minimal() +
  geom_vline(xintercept=c(1.5, 2.5, 3.5, 4.5, 5.5)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_text(data=counts, aes(x=as.factor(slack_prop), 
                             label=scales::percent(median_fails, accuracy=0.1), 
                             colour=as.factor(sample_sites), ymax=3, y = q90_y_pos), 
            angle=-90, 
            position = position_dodge(width = 1), show.legend = FALSE)

