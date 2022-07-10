#-----------------------------------------------------------------------#
#                               Analysis script
# Temporal Dynamics of depressive Symptomatology: 
# An idiographic time series analysis
# 
# Created by: Björn Siepe
# Last edited: 30.11.2021
#                                   Part 5
#                               Block Bootstrap
#-----------------------------------------------------------------------#

# Description -------------------------------------------------------------
# This script is used for running the block bootstrap on our server

# Packages and data -------------------------------------------------------
if(!require(mgm)){
  install.packages("mgm", dependencies = TRUE)
}
library(mgm)
# Set wd
setwd("/mnt/md0/siepe/")

# Load tvvar results
load("tvvar_res.RData")
load("data.RData")

# Create id vector manually
pps <- c(1,2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21,22)
daily_vars <- c("phq1", "phq2", "sleep", 
                "rumin", "soc_quant", "soc_qual") 

# Block Bootstrap ---------------------------------------------------------
# Performing a non-parametric block bootstrap for each estimate
# Set up list to store results in
rel_res <- list()

# Start bootstrap
set.seed(2021)
for(i in pps){
  data_ind_rel <- subset(data, id == i)                      # same as above: get individual data
  data_ind_rel <- data_ind_rel[daily_vars]                   # only use self-report variables
  rel_res[[i]] <- mgm::resample(object = tvvar_res[[i]],         # choose individual results that were obtained above for resampling
                            data = data_ind_rel, 
                            nB = 1000,                             # 1000 bootstrap samples
                            blocks = 20,                           # 20 blocks for nonparametric block bootstrap
                            seeds = 1:1000,                        # seeds used to sample 
                            quantiles = c(.05, .95))               # obtain 5% and 95% quantiles of bootstrapped sampling distribution
  rel_res[[i]]$models <- NULL                       
  }

# Output ------------------------------------------------------------------
save(rel_res, file = "Bootstrapping Results.Rdata")
saveRDS(rel_res, file = "Bootstrapping Results2.rds")


