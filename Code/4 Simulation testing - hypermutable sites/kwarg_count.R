#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))

dat <- read.table(args[1], header=TRUE, fill=TRUE) %>%
  filter(!is.na(R)) %>%
  filter(R == 0) %>%
  mutate(Pars = RM + SE)
cat(paste0(" ", min(dat$Pars), "\n"))
