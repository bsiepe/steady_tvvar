#-----------------------------------------------------------------------#
#                               Analysis script
# Temporal Dynamics of depressive Symptomatology: 
# An idiographic time series analysis
# 
# Created by: Björn Siepe
#                                   Part 3
#                               Main Analysis
#-----------------------------------------------------------------------#

# Description -------------------------------------------------------------
# In this script, the following steps are performed: 
# 1. Bandwidth Selection + Analysis
# 2. Model estimation with time-varying VAR 
# 3. Obtaining results of bootstrapped sampling distribution from server script
# 4. Computing Predictability
# 5. Visualization
  # 5.1. Network visualizations
  # 5.2. Time-Varying parameters
# 6. Test for stationarity
# 7. Network plot animation

# NOTES #
# This code builds heavily on the very helpful tutorial paper + Code
# by Haslbeck & Bringmann (2021)
# https://github.com/jmbh/tvvar_paper/blob/master/Tutorials/tutorial_mgm.R


# -------------------------------------------------------------------------
# Loading data/packages & preparations ------------------------------------
# -------------------------------------------------------------------------
library(mgm)          # modeling
library(qgraph)       # network plots
library(ggplot2)      # plotting
library(dplyr)        # data manipulation
library(mlVAR)        # data simulation
library(gridExtra)    # arranging plots in grid
library(flextable)    # creating nice APA tables
library(here)         # makes sourcing easier
library(stats)        # calculating median

# Create a data directory
data_dir <- "Data"

# Create a figure directory
figure_dir <-"Output/"

# Load  data
load(paste0(data_dir, "/Processed/data.RData"))

# Set up the loops by creating a vector with unique IDs 
pps <- unique(data$id)   

# -------------------------------------------------------------------------
# Bandwidth selection -----------------------------------------------------
# -------------------------------------------------------------------------
bw_res <- list()         # store full results of bandwidth selection

# Specify the bandwidth selection sequence to choose from
bw_seq <- seq(0.01, 1, length = 10)

# Create list of relevant variables
# Note: This changes the column order in data
# This is the column order that will be used for the rest of the analyses
daily_vars <- c("phq1", "phq2", "sleep", 
                "rumin", "soc_quant", "soc_qual") 

# Perform bandwidth selection for each ID 
# estimate a model at each point in the test set, weighing the data at the 
# points in the test set with 0
set.seed(2021)
time_before_bw <- Sys.time()               # Check how long it takes
for(i in pps){                 
  data_ind_bw <- subset(data, id == i)     # only use data with current ID
  data_ind_bw <- data_ind_bw[daily_vars]   # relevant variables for the model
  bw_res[[i]] <- mgm::bwSelect(data = data_ind_bw,
                               type = rep("g", 6), # 6 Gaussian variables
                               level = rep(1, 6),  # 6 continuous variables (indicated by level = 1)
                               bwSeq = bw_seq,     # sequence of bandwidths
                               bwFolds = 5,        # 5-fold CV
                               bwFoldsize = nrow(data_ind_bw)/5, # equally sized folds
                               modeltype = "mvar",
                               lags = 1,           # use lag 1
                               scale = TRUE,       # scale the variables
                               pbar = TRUE)        # use progress bar
  
}
time_after_bw <- Sys.time()-time_before_bw 
# Or load results directly: 
# load(file = paste0(data_dir, "/Processed/bw_results.Rdata"))

# -------------------------------------------------------------------------
# Optimal bandwidth for every id ------------------------------------------
# -------------------------------------------------------------------------
# Goal: Obtain bandwidth that minimizes the mean error over the test set (called bw)
# and the associated error (called mse)
# Set up list and loop through results
min_bw <- list()
for(i in pps){
  min_bw[[i]] <- data.frame(id = i,        
                            bw = bw_seq[which.min(bw_res[[i]]$meanError)],
                            mse = min(bw_res[[i]]$meanError))
}
# Now combine the results to get a dataframe with id and minimal bandwidth
df_min_bw <- do.call(rbind, min_bw)

# Second BW selection for extreme cases -----------------------------------
# For min_bw = 0.01 or = 1, we do another round of bw selection
# Note: the largest value for bw_seq_small is 0.009, not 0.0091, due to rounding
# should not make much of a difference and was chosen arbitrarily anyway
bw_seq_small <- round(seq(0.001, 0.0091, length = 10),3)
bw_seq_large <- round(seq(1.1, 2, length = 10),3) 

# Find individuals with bw = 0.01
# extract their ids for looping
bw_small_ind <- df_min_bw %>% 
  filter(bw == 0.01) %>%
  pull(id)

# Find individuals with bw = 1 
# extract their ids for looping
bw_large_ind <- df_min_bw %>% 
  filter(bw == 1) %>%
  pull(id)

### Issue: id = 10 and id 19 will throw an error in the second iteration:
# "Error in elnet(xd, is.sparse, ix, jx, y, weights, offset, type.gaussian,  : 
# y is constant; gaussian glmnet fails at standardization step"
# This is likely due to floor/ceiling effects wihtin either phq1 or soc_qual
data %>% 
  filter(id == 10) %>% 
  ggplot(aes(x = seq(1,316, by = 1), y = phq1))+geom_line()
# These are problematic when testing small bandwidths due to small windows that
# then have constant response data 

# Code to test this:
# set.seed(2021)
# for(i in 10){
#   data_ind_bw <- subset(data, id == i)     # only use data with current ID
#   data_ind_bw <- data_ind_bw[daily_vars]   # relevant variables for the model
#   bw_res_small[[i]] <- mgm::bwSelect(data = data_ind_bw,
#                                      type = rep("g", 6), # 6 Gaussian variables
#                                      level = rep(1, 6),  # 6 continuous variables
#                                      bwSeq = bw_seq_small,   # use new small bw_seq
#                                      bwFolds = 5,
#                                      bwFoldsize = nrow(data_ind_bw)/5, 
#                                      modeltype = "mvar",
#                                      lags = 1,     # use lag 1
#                                      scale = TRUE,  # scale the variables
#                                      pbar = TRUE) # use progress bar
#   
# } 

# Therefore, exclude id = 10 and 19 for now
bw_small_ind <- bw_small_ind[!bw_small_ind %in% c(10,19)]

### New Loop for small bandwidths
bw_res_small <- list()

# Start the loop with the smaller bw sequence
set.seed(2021)
for(i in bw_small_ind){
  data_ind_bw <- subset(data, id == i)     # only use data with current ID
  data_ind_bw <- data_ind_bw[daily_vars]   # relevant variables for the model
  bw_res_small[[i]] <- mgm::bwSelect(data = data_ind_bw,
                               type = rep("g", 6), # 6 Gaussian variables
                               level = rep(1, 6),  # 6 continuous variables
                               bwSeq = bw_seq_small,   # use new small bw_seq
                               bwFolds = 5,
                               bwFoldsize = nrow(data_ind_bw)/5, 
                               modeltype = "mvar",
                               lags = 1,     # use lag 1
                               scale = TRUE,  # scale the variables
                               pbar = TRUE) # use progress bar
  
} 

# Or load results directly: 
# load(file = paste0(data_dir, "/Processed/bw_results_small.RData"))

### Do the same for large bandwidths
bw_res_large <- list()
set.seed(2021)
for(i in bw_large_ind){
  data_ind_bw <- subset(data, id == i)     # only use data with current ID
  data_ind_bw <- data_ind_bw[daily_vars]   # relevant variables for the model
  bw_res_large[[i]] <- mgm::bwSelect(data = data_ind_bw,
                                     type = rep("g", 6), # 6 Gaussian variables
                                     level = rep(1, 6),  # 6 continuous variables
                                     bwSeq = bw_seq_large,   # use new large bw_seq
                                     bwFolds = 5,
                                     bwFoldsize = nrow(data_ind_bw)/5, 
                                     modeltype = "mvar",
                                     lags = 1,     # use lag 1
                                     scale = TRUE,  # scale the variables
                                     pbar = TRUE) # use progress bar
  
} 
# Or load results directly
# load(file = paste0(data_dir, "/Processed/bw_results_large.RData"))

## Check if there was an improvement in bandwidth
# As above, obtain bandwidths + their respective errors
min_bw_small <- list()
for(i in bw_small_ind){
  min_bw_small[[i]] <- data.frame(id = i,        
                            bw = bw_seq_small[which.min(bw_res_small[[i]]$meanError)],
                            mse = min(bw_res_small[[i]]$meanError))
}
# get dataframe
df_min_bw_small <- do.call(rbind, min_bw_small)

## same for large BWs
min_bw_large <- list()
for(i in bw_large_ind){
  min_bw_large[[i]] <- data.frame(id = i,        
                                  bw = bw_seq_large[which.min(bw_res_large[[i]]$meanError)],
                                  mse = min(bw_res_large[[i]]$meanError))
}
df_min_bw_large <- do.call(rbind, min_bw_large)

# attach them both together for all IDs with a second iteration
df_min_bw_second <- rbind(df_min_bw_small, df_min_bw_large)

## Check if MSE of new bandwidth is smaller
# If MSE of new bandwidth is smaller than MSE of original bandwidth, 
# choose the new one
# Final bandwidth is saved under bw.final
df_min_bw <- df_min_bw %>% 
  left_join(df_min_bw_second, by = "id") %>% 
  mutate(bw_final = if_else(mse.y >= mse.x | is.na(mse.y), bw.x, bw.y)) %>% 
  mutate(mse_final = if_else(mse.y >= mse.x | is.na(mse.y), mse.x, mse.y)) %>% 
  select(id, bw_final, mse_final)

# Analyze bandwidth selection ---------------------------------------------
# Inspect the bandwidths and respective error for each id
# this only works for IDs with the original bw_seq
for(i in pps){
  plot(x = bw_res[[i]]$call$bwSeq,
       y = bw_res[[i]]$meanError,
       main = paste0("Participant ", i),
       xlab = "Bandwidth",
       ylab = "MSE",
       ylim = c(min(bw_res[[i]]$meanError)-0.05,
                max(bw_res[[i]]$meanError)+0.05),
       # Add numerical errors to each point
       text(bw_res[[i]]$call$bwSeq,
            y = bw_res[[i]]$meanError,
            labels = round(bw_res[[i]]$meanError,4),
            cex = 0.5,
            pos = 3))
  # Plot a line for the minimum error
  abline(h = min(bw_res[[i]]$meanError), lty = 2)
}

### Compute average effectively used sample size 
# Loop over the number of effective sample size for each timepoint for each id
# For example, a time series length of 100 and Ne of 50 means
# that half of the information in the time series was used
m_eff_ss <- list()
for(i in pps){
  m_eff_ss[[i]] <- data.frame(id = i, mean_ss = mean(tvvar_res[[i]]$Ne))
}
m_eff_ss <- do.call(rbind, m_eff_ss)

# Attach individual bandwidth
m_eff_ss %>% 
  left_join(df_min_bw, by = "id")

# Illustrate different bandwidths --------------------------------------
# For Figure 1 in the manuscript, illustrate the difference
# in weighting of timepoints of a small, medium and a larger bandwidth
# Use timepoint = 10 for plotting (roughly middle of the time series)

# Save weights of participant 2 (small bw)
df_bw_weights_02 <- data.frame(
  tp = seq(1,539, by = 1),
  weights = tvvar_res[[2]]$tvmodels[[10]]$call$weights) 
# Save weights of participant 3 (medium bw)
df_bw_weights_03 <- data.frame(
  tp = seq(1,308, by = 1),
  weights = tvvar_res[[3]]$tvmodels[[10]]$call$weights)
# Save weights of participant 21 (large bw) 
df_bw_weights_21 <- data.frame(
  tp = seq(1,194, by = 1),
  weights = tvvar_res[[21]]$tvmodels[[10]]$call$weights)

# Plot as line plots in a grid
# First plot
tmp_bw_plot <- gridExtra::grid.arrange(ggplot(df_bw_weights_02, aes(x = tp, y = weights))+
  geom_point(size = 0.3,alpha = 0.5)+
  # area under the curve
  geom_ribbon(aes(ymin = 0, ymax =  weights), fill = 'lightblue')+
  # area above the curve
  geom_ribbon(aes(ymin = weights, ymax = 1), fill = "grey", alpha = 0.4)+
  
  # Themes/Colours
  theme_bw()+
  theme(legend.text = element_text(size = 8),
        plot.title = element_text(size = 8),
        panel.grid.major.x = element_blank(), # remove gridlines
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", size = 1))+
  guides(fill = "none")+
  # Labels
  labs(x = "Day", y  = "Weight",
       title = "ID 2, BW = 0.009")+
  # Scales
  scale_x_continuous(breaks = c(1,100,200,300,400,500),
                     expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)),
  
  # second plot
  ggplot(df_bw_weights_03, aes(x = tp, y = weights))+
    geom_point(size = 0.3,alpha = 0.5)+
    # area under the curve
    geom_ribbon(aes(ymin = 0, ymax =  weights), fill = 'lightblue')+
    # area above the curve
    geom_ribbon(aes(ymin = weights, ymax = 1), fill = "grey", alpha = 0.4)+
    
    # Themes/Colours
    theme_bw()+
    theme(legend.text = element_text(size = 8),
          plot.title = element_text(size = 8),
          panel.grid.major.x = element_blank(), # remove gridlines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", size = 1))+
    guides(fill = "none")+
    # Labels
    labs(x = "Day", y  = "Weight",
         title = "ID 3, BW = 0.12")+
    # Scales - 
    scale_x_continuous(breaks = c(1,100,200,300),
                       expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)),
  
  # third plot
  ggplot(df_bw_weights_21, aes(x = tp, y = weights))+
    geom_point(size = 0.3,alpha = 0.5)+
    # area under the curve
    geom_ribbon(aes(ymin = 0, ymax =  weights), fill = 'lightblue')+
    # area above the curve
    geom_ribbon(aes(ymin = weights, ymax = 1), fill = "grey", alpha = 0.4)+
    
    # Themes/Colours
    theme_bw()+
    theme(legend.text = element_text(size = 8),
          plot.title = element_text(size = 8),
          panel.grid.major.x = element_blank(), # remove gridlines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", size = 1))+
    guides(fill = "none")+
    # Labels
    labs(x = "Day", y  = "Weight", 
         title = "ID 21, BW = 1.00")+
    # Scales 
    scale_x_continuous(breaks = c(1,100,200),
                       expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)),
  nrow = 1)

ggsave(tmp_bw_plot, filename = paste0(figure_dir,"bw_comparison.svg"),
       device = "svg", width = 16, height = 9, units = "cm")



# -------------------------------------------------------------------------
# Model estimation --------------------------------------------------------
# -------------------------------------------------------------------------
# Build time-varying VAR model for every ID with bandwidths obtained above

# Set up an empty list to store results
tvvar_res <- list()  

# Loop over all ids
time_before_model <- Sys.time()

# Set seed again and start modelling
set.seed(2021)
for(i in pps){
  data_ind <- subset(data, id == i)                               # create individual data for i
  data_ind <- data_ind[daily_vars]                                # only select self-report variables
  tvvar_res[[i]] <- mgm::tvmvar(data = data_ind,                         
                                type = rep("g", 6),                    # all gaussian variables
                                level = rep(1, 6),                     # all continuous variables
                                lambdaSel = "EBIC",                    # use EBIC for selection of regularization parameter
                                lambdaGam = 0,                         # set EBIC hyperparameter
                                lags = 1,                              # choose VAR(1) model
                                consec = 1:nrow(data_ind),             # indicate which timepoints are consecutive
                                scale = TRUE,                          # scale all variables (important bc of regularization)
                                estpoints = seq(0, 1, length = 20),    # estimation points
                                bandwidth = df_min_bw[df_min_bw$id == i,]$bw_final, # use bandwidth that was selected before for each id
                                pbar = TRUE)                           # progress bar to show if it is worth to get another coffee
  
}
# Check how long this took
time_after_model <- Sys.time()-time_before_model

# Or load results directly:
# load(file = paste0(data_dir, "/Processed/tvvar_res.RData"))


### Compute average effectively used sample size 
# Loop over the number of effective sample size for each timepoint for each id
# For example, a time series length of 100 and Ne of 50 means
# that half of the information in the time series was used
m_eff_ss <- list()
for(i in pps){
  m_eff_ss[[i]] <- data.frame(id = i, mean_ss = mean(tvvar_res[[i]]$Ne))
}
m_eff_ss <- do.call(rbind, m_eff_ss)

# Attach individual bandwidth
m_eff_ss %>% 
  left_join(df_min_bw, by = "id")

### Inspect effectively used sample size for IDs 2,3,21 and estimation point 10
tvvar_res[[2]]$Ne[10]
tvvar_res[[3]]$Ne[10]
tvvar_res[[21]]$Ne[10]

# -------------------------------------------------------------------------
# Bootstrapped sampling distribution --------------------------------------
# -------------------------------------------------------------------------
# Perform a non-parametric block bootstrap to obtain bootstrapped
# sampling distributions of parameters

# Empty list to store the results is
rel_res <- list()                                                

# WARNING: Doing this for the whole sample takes ~ 10 days on my computer
# Save the time that this takes
time_before_rel <- Sys.time()
# Compute the reliability of each estimate for every individual
# Performing a non-parametric block bootstrap
for(i in pps){
  data_ind_rel <- subset(data, id == i)                      # same as above: get individual data
  data_ind_rel <- data_ind_rel[daily_vars]                   # only use self-report variables
  rel_res[[i]] <- mgm::resample(object = tvvar_res[[i]],         # choose individual results that were obtained above for resampling
                            data = data_ind_rel, 
                            nB = 1000,                             # 1000 bootstrap samples
                            blocks = 20,                           # 20 blocks for nonparametric block bootstrap
                            seeds = 1:1000,                        # seeds used to sample
                            quantiles = c(.05, .95))               # obtain 5% and 95% quantiles of bootstrapped sampling distribution
  rel_res[[i]]$models <- NULL     # delete the 1000 model objects, otherwise this takes up way too much disk space
}

# Compute how long it toook
time_after_rel <- Sys.time()-time_before_rel


# Reading in all bootstrap server results ---------------------------------
# rel_res <- readRDS(file = paste0(data_dir, "/Processed/BootstrappingResults2.rds"))


# -------------------------------------------------------------------------
# Predictability ----------------------------------------------------------
# -------------------------------------------------------------------------
# Compute the Predictability of each node
# again, start with an empty list
pred_res <- list()

# Compute R-Squared and Root Mean Squared Error for each node for each ID
for(i in pps){
  data_ind_pred <- subset(data, id == i)                # use data for one id
  data_ind_pred <- data_ind_pred[daily_vars]            # only use relevant variables
  pred_res[[i]] <- predict(object = tvvar_res[[i]],    
                           data = data_ind_pred, 
                           errorCon = c("R2", "RMSE"),      # relevant metrics
                           tvMethod = "weighted",           # weighted average of predictions of all models
                           consec = 1:nrow(data_ind_pred))  # important: also update this when consec does not work
  
}

## Save average R2 and RMSE for each individual
# Iterate over list of prediction results, save in another list
l_pred_errors <- list()
for(i in pps){
  l_pred_errors[[i]]<- pred_res[[i]]$errors
  l_pred_errors[[i]]$id <- i
}
# Save as dataframe
df_pred_errors <- do.call(rbind, l_pred_errors)
# Mean and sd of R-Squared across IDs and variables
mean(df_pred_errors$R2)
sd(df_pred_errors$R2)
# Mean and sd of RMSE across IDs and variables
mean(df_pred_errors$RMSE)
sd(df_pred_errors$RMSE)
# Obtain errors per ID
df_pred_errors %>% 
  group_by(id) %>% 
  summarize(m_RMSE = mean(RMSE),
            m_R2 = mean(R2))

# Small visualization of R-Squared averaged over all vars
library(ggbeeswarm)
ggplot(df_pred_errors, aes(color = Variable))+
  geom_boxplot(aes(x = Variable, y = R2))+
  coord_flip()+
  geom_beeswarm(aes(x = Variable, y = R2, color = Variable), cex = 2.5)+
  theme_classic()

# -------------------------------------------------------------------------
# Visualization -----------------------------------------------------------
# -------------------------------------------------------------------------
### Goals:
# 1. Get Model visualization at timepoints 2, 10, 19
# 2. Plot specific parameters over time, find those with highest changes
# 3. Store these graphics neatly

# First, get out the parameters and get the correct sign for them 
# ind_negative serves as an indicator for every parameter that has a negative sign
tvvar_pars <- list()
ind_negative <- list()
for(i in pps){
  # get weights
  tvvar_pars[[i]] <- tvvar_res[[i]]$wadj
  # indicator for parameters that have a negative sign, gives position in array
  ind_negative[[i]] <- which(tvvar_res[[i]]$signs == -1, arr.ind = T)
  # give parameters their correct negative sign
  tvvar_pars[[i]][ind_negative[[i]]] <- tvvar_pars[[i]][ind_negative[[i]]] * -1 
}


# Network plots -----------------------------------------------------------
### Plot network at timepoint 2, 10, 19 (beginning, middle, end)
# Do this for every individual and save the results as graphic

# Name the variables for network plot
var_labels_n <- c("Anhe-\ndonia", "Feeling \ndown", "Sleep",
                  "Rumi-\nnation", "Social \nQuantity", "Social \nQuality")


## The following for loop creates 3 network plots for each 
# and stores them in a single pdf  
# Do not plot the outermost ones(1,20) as these might use too few data
# Transposing the weight matrix needs to be done as qgraph and mgm
# have different ways to store edges (see Haslbeck & Waldorp, 2020, p.30)
for(i in pps){
  svg(filename = paste0(figure_dir,"3network_id", i,".svg"),
      width = 16, height = 9)
  par(mfrow=c(1,3))
  for(ep in c(2, 10, 19)){                                       # choose 3 timepoints
    print(qgraph(t(tvvar_pars[[i]][, , 1, ep]),       # transpose the weights matrix for qgraph default 
                 layout = "circle",                                    # choose layout of network
                 labels = var_labels_n,                                # proper labelling
                 edge.color = t(tvvar_res[[i]]$edgecolor[, , 1, ep]),  # transpose for qgraph and get edgecolor out of model object
                 mar = rep(5, 4),                                      # margins of the graph
                 vsize=15, esize=7, asize=6,                           # adjust variable size
                 maximum = .5,                                         # max. effect for plotting relations
                 pie = pred_res[[i]]$tverrors[[ep]][, 3],              # add prediction error to each node
                 title = ep,                                           # est.point as title
                 title.cex = 5,                                        # large title font
                 negDashed = TRUE,
                 label.fill.horizontal = 0.9))                         # font filling respective to node size
    
    box("figure", bty = "o", col = "grey") # add box around plots
    }
 dev.off()
}


### Save all networks at all timepoints in one pdf
pdf(file = paste0(data_dir, "/Processed/network_plots_all_ids.pdf"))
for(i in pps){
  par(mfrow=c(7,3)) 
  for(ep in c(1:20)){                                       # choose 3 timepoints
    print(qgraph(t(tvvar_pars[[i]][, , 1, ep]),       # transpose the weights matrix for qgraph default 
                 layout = "circle",                                    # choose layout of network
                 labels = var_labels_n,                                # proper labelling
                 edge.color = t(tvvar_res[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
                 mar = rep(5, 4),                                      # margins of the graph
                 vsize=20, esize=15, asize=13,                         # adjust variable size
                 maximum = 0.5,                                        # max. effect for plotting relations 
                 pie = pred_res[[i]]$tverrors[[ep]][, 3],              # add prediction error to each node
                 title = paste0("Estimation point = ", ep, ", ID = ", i),        
                 title.cex=0.9,
                 label.fill.horizontal = 0.9))
    
    box("figure", bty = "o", col = "grey") # add box around plots
    }
}
dev.off()

### Investigate 0.5 as maximum for edge weight plotting
# Check how many absolute non-zero edge weights were below 0.5
nonzero_ests <- list()
for(i in pps){
  nonzero_ests[[i]] <- tvvar_res[[i]]$wadj[tvvar_res[[i]]$wadj != 0]
}
m_nonzero_ests <- do.call(rbind, nonzero_ests)
sum(m_nonzero_ests <= 0.5)/length(m_nonzero_ests)


# Network plot ID 2 -------------------------------------------------------
# For ID 2, create a plot with all non-empty networks to illustrate extent of variation
# over time
for(i in 2){
  svg(filename = paste0(figure_dir,"ppt_networks_id", i,".svg"),
      width = 16, height = 18)
  par(mfrow=c(2,3)) 
  for(ep in c(1,2,3,12,14,15)){                           # choose non-empty estimationpoints
    print(qgraph(t(tvvar_pars[[i]][, , 1, ep]),       # transpose the weights matrix for qgraph default 
                 layout = "circle",                                    # choose layout of network
                 labels = var_labels_n,                                # proper labelling
                 edge.color = t(tvvar_res[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
                 mar = rep(5, 4),                                      # margins of the graph
                 vsize=15, esize=7, asize=6,                           # adjust variable size
                 maximum = 0.5,                                        # max. effect for plotting relations 
                 pie = pred_res[[i]]$tverrors[[ep]][, 3],              # add prediction error to each node
                 title = ep,                                           # estimation point as title        
                 title.cex=5,                                          # title size
                 negDashed = TRUE,
                 label.fill.horizontal = 0.9))                         # label size
    
    box("figure", bty = "o", col = "grey") # add box around plot
    }
}
dev.off()

### When did significant life events happen?
# As per the file "Life_Events.xlsx", two important events:
# 7/18: increased stress at work for the whole month
# 10/18: severe medical diagnosis of a friend
# As we do not know when the diagnosis happened during the month, 
# we take the 15th of october as the middle of the month
# IMPORTANT TODO: These specific dates will need to be obscured when publishing 
# the Code online

# Obtain relative position of events in the time series
tmp_data_02 <- data %>% filter(id == 2)
# Divide row number of dates by total number of observations
which(tmp_data_02$date == "2018-07-01")/nrow(tmp_data_02)
which(tmp_data_02$date == "2018-10-15")/nrow(tmp_data_02)
# Compare with estimation points (also denoted relative to full length of time series)
tvvar_res[[2]]$call$estpointsNorm

### Average edge strength of present edges
# Find all edges that were non-zero and look at their size
# use median because of 1 large outlier
stats::median(tvvar_res[[2]]$wadj[tvvar_res[[2]]$wadj != 0])


# Plot parameters over time -----------------------------------------------
### GOAL: Plot time-varying parameters over time
# Again, especially this part of the code builds heavily on Haslbeck & Bringmann (2021)
# It is simply adapted to a n>1 setting

## Goal: Compute the standard deviation for every parameter for every individual
# First calculate SDs
# Then properly transform the matrices

# set up a loop to get a 6x6 matrix of SDs for each individuals
pars_sds <- list()
for(i in pps){
  pars_sds[[i]] <- apply(tvvar_pars[[i]], c(1,2), sd)   # second argument: apply this to rows and columns
}

## Goal: matrix for each individual that has pxp rows and three columns
# first column: variable that is affected (dv)
# second column: variable from which the effect comes (iv)
# third column: standard deviation

# First set up a loop to get an empty matrix for each individual
# 6^2 for results of 6 parameters
pars_sds_mat <- list()
for(i in pps){
  pars_sds_mat[[i]] <- matrix(NA, 6^2, 3)
}

# Loop to fill the matrix with the standard deviations
# k is an index for the row number
# j is an index for the column number
for(i in pps){
  counter <- 1
  for(k in 1:6){
    for(j in 1:6) {
      pars_sds_mat[[i]][counter, ] <- c(k, j, pars_sds[[i]][k, j]) 
      counter <- counter + 1
    }
  }
}
# pars_sds_mat now contains the standard deviations for each parameter

### Goal: Select 3 parameters with highest SDs for each individual for plotting
# To do so, first order pars_sds_mat to obtain highest SDs
# Convert to dataframes here, makes it more convenient (esp. for later plotting looping)
pars_sds_mat_df <- lapply(pars_sds_mat, as.data.frame)

# Order data frames by SD
# need is.vector() as some list entries are NULL (due to IDs being excluded)
pars_sds_mat_df <- lapply(pars_sds_mat_df, function(x) if(is.vector(x$V3)) x[order(x$V3, decreasing = TRUE), ] else NULL) 

# Only get 3 parameters with highest SDs
pars_sds_mat_df <- lapply(pars_sds_mat_df, function(x) head(x, 3))

## Goal: make line plots with uncertainty bands
# Also, automatize naming of the parameters in the plot
# First, make a matrix to indicate which parameters will be used and store in list
# call this list pars_disp
# V1 is the DV, V2 is the IV   
pars_disp <- list()
for(i in pps){
  pars_disp[[i]] <- as.matrix(pars_sds_mat_df[[i]][1:3,1:2], nrow = 2, byrow = TRUE)
}

# Save proper variable descriptions for plotting
# This needs to be double checked! (right order or not)
plot_vars <- c("Anhedonia", "Feeling Down","Sleep", "Rumination",
               "Social Quant.", "Social Qual." ) 

####### BEGINNING Loop
# Get a data frame with pointests and CIs for every parameter
# repeat the variable number 20 times each (one row for each estimation point)
pars_plot <- list()
for(i in pps){
  pars_plot[[i]] <- data.frame(dv = rep(pars_disp[[i]][,1], each = 20),
                               iv = rep(pars_disp[[i]][,2], each = 20),
                               ep = seq(1,20,1),    # get estimationpoint from 1 to 20 
                               pointest = rep(NA, 60),
                               ci_low = rep(NA, 60),
                               ci_up = rep(NA, 60))
}

# Now extract cis and point estimates for every person through a loop
# pars_plot then contains relevant parameters at every timepoint with cis and estimates
# Uncomment CI part when sampling distribution should be visualized
# again, these are not true CIs around a parameter 
low <- list()
up <- list()
pointest <- list()
for(i in pps){
  for(j in 1:3){  # 3 parameters
    # # anonymous function to get out the 5% and 95% quantile
    # cis <- apply(rel_res[[i]]$bootParameters[pars_disp[[i]][j,1], pars_disp[[i]][j,2], 1, , ], 1,
    # function(x) {quantile(x, probs = c(.05, .95))})
    # # # save upper and lower boundary
    # low[[j]] <- cis[1,]
    # up[[j]] <- cis[2,]
    # get the point estimate and store it
    est <- tvvar_pars[[i]][pars_disp[[i]][j,1], pars_disp[[i]][j,2], 1, ]
    pointest[[j]] <- est
  }
  pars_plot[[i]]$pointest <- unlist(pointest)
   # pars_plot[[i]]$ci_low <- unlist(low)
   # pars_plot[[i]]$ci_up <- unlist(up)
}

# Now transform these data frames so that dv and iv come into one column
# name it as "iv on dv", so for example "4 on 2" 
# for an effect from parameter 4 to parameter 2
# This is helpful for grouping/plotting
pars_plot <- lapply(pars_plot, function(x){
  x$effect <- as.factor(paste(x$iv, x$dv, sep=" on ")) # concatenate the columns  
  x              # return data.frame for proper storing
})

### Build the Plot
# Choose colors manually
# Took these from viridis plasma colour palette
# see https://waldyrious.net/viridis-palette-generator/
plot_colors <- c('#f0f921',  '#cc4778', '#0d0887')

## Set up the loop for ggplot
# Uncomment geom_ribbon for a display of bootstrap results
# This sometimes means that ylim needs to be adapted for large sampling dists.
for(i in pps){
  # get names of variables for legend, store them temporarily
  dv1 <- plot_vars[pars_disp[[i]][1,1]]
  iv1 <- plot_vars[pars_disp[[i]][1,2]]
  dv2 <- plot_vars[pars_disp[[i]][2,1]]
  iv2 <- plot_vars[pars_disp[[i]][2,2]]
  dv3 <- plot_vars[pars_disp[[i]][3,1]]
  iv3 <- plot_vars[pars_disp[[i]][3,2]]
  # Start plotting
  tmp_plot <- ggplot(data = pars_plot[[i]], aes(x = ep, 
                                                colour = effect, fill = effect))+
    # Setting geoms
    geom_line(aes(y = pointest, linetype = effect), size = 1)+
    # geom_ribbon(aes(ymin=ci_low, ymax = ci_up), alpha = 0.3, linetype = 0)+
    
    # Themes/Colours
    theme_light()+
    theme(legend.text = element_text(size = 8, margin = margin(r = 25, unit = "pt")), #space between text and next key
          legend.position="bottom",          # horizontal legend
          legend.spacing.x = unit(4, 'pt'),  # spacing between legend keys
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,0,-10), # bring legend closer to figure
          panel.grid.major.x = element_blank(), # remove gridlines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),       # remove borders
          axis.line = element_line(colour = "black"))+  # use axis lines
    
    # now set all colors and labels manually
    # the bquote expression uses the variables stored above for the legend
    # see https://lukemiller.org/index.php/2017/05/r-plotmath-functions-combined-with-variable-values/
    scale_colour_manual(values = plot_colors,
                        labels = c(bquote(.(iv1)["t-1"]*''%->%''*.(dv1)["t"]), 
                                   bquote(.(iv2)["t-1"]*''%->%''*.(dv2)["t"]), 
                                   bquote(.(iv3)["t-1"]*''%->%''*.(dv3)["t"])))+
    scale_fill_manual(values = plot_colors,
                      labels = c(bquote(.(iv1)["t-1"]*''%->%''*.(dv1)["t"]), 
                                 bquote(.(iv2)["t-1"]*''%->%''*.(dv2)["t"]), 
                                 bquote(.(iv3)["t-1"]*''%->%''*.(dv3)["t"])))+
    scale_linetype_manual(values = c(1,2,3),
                          labels = c(bquote(.(iv1)["t-1"]*''%->%''*.(dv1)["t"]), 
                                     bquote(.(iv2)["t-1"]*''%->%''*.(dv2)["t"]), 
                                     bquote(.(iv3)["t-1"]*''%->%''*.(dv3)["t"])))+
    
    # Axes
    scale_x_continuous(breaks = c(1,2,5,10,15,19,20), 
                       limits = c(1,20),
                       expand = c(0,0))+   # makes plot end with axis
    scale_y_continuous(breaks = c(-0.50, -0.25, 0.00, 0.25, 0.50),
                       limits = c(-0.5,0.5),
                       expand = c(0,0))+   # makes plot end with axis
    
    # Description
    labs(x = "Estimation Point",
         y = "Estimate",
         colour = "Effect",
         fill = "Effect",
         linetype = "Effect")
  ggsave(tmp_plot, filename = paste0(figure_dir,"varying_pars_id_", i, ".svg"),
         device = "svg", width = 16, height = 9, units = "cm")
  
}

# Time-Varying Parameter plot ID 6 ----------------------------------------
# For Participant 6, we want to plot 1 additional variable
# and annotate some life events
# This code is not made available publicly because it hardcodes
# the dates of certain critical life events. The code can be shared
# upon request


# -------------------------------------------------------------------------
# Test for time-varying ---------------------------------------------------
# -------------------------------------------------------------------------
# This test is described in the manuscript
# In short: We want to compare the empirical error of the actual model
# against a sampling distribution
# Sampling distribution is obtained by fitting a normal VAR model 
# and then simulating 100 models with the obtained parameters
# On each of these, we fit a tvvar model and compute the RMSE
# These 100 RMSEs then serve as sampling distribution

# Deviations from the computation in the tutorial by Haslbeck:
# 1. Use EBIC instead of CV for var and tv-var models
# 2. Use consec instead of timepoints in estimation of tv-var because we have no missings
# 3. Use tvmethod = "weighted" for prediction errors (as we use above)

# Get RMSE from time-varying model ----------------------------------------
error_emp <- list()
for(i in pps){
  error_emp[[i]] <- mean(pred_res[[i]]$errors$RMSE)    
}

# Fit stationary model ----------------------------------------------------
# Here, we fit a simple stationary mvar model on the timeseries of each individual
var_res <- list()
set.seed(2021)
for(i in pps){
  data_ind <- subset(data, id == i) # create individual data for i
  data_ind <- data_ind[daily_vars]  # use only relevant variables
  var_res[[i]] <- mgm::mvar(data = data_ind,
                            type = rep("g", 6),
                            level = rep(1, 6), 
                            lambdaSel = "EBIC",      
                            lambdaGam = 0,
                            lags = 1,
                            consec = 1:nrow(data_ind),
                            pbar = TRUE)
  
}

# Simulate data and fit time-varying model to it --------------------------
# This is a way to extract the no. of included data points
# for each individual, which is relevant for the simulation
sum(tvvar_res[[2]]$call$data_lagged$included) 

nIter <- 100 # number of iterations to draw distribution
l_data <- list()   # set up empty lists as storage
l_model <- list()
l_meanerror <- rep(NA, nIter)

# get parameters out of each model object
pars_sim <- list()
intercepts_sim <- list()

# loop to get parameters out
for(i in pps){
  # get out parameter weights ([,,1] means lag-1)
  pars_sim[[i]] <- var_res[[i]]$wadj[,,1]
  # give negative parameters a negative sign
  pars_sim[[i]][var_res[[i]]$edgecolor=="red"] <- pars_sim[[i]][var_res[[i]]$edgecolor=="red"]*-1
  # also get out the intercepts
  intercepts_sim[[i]] <- unlist(var_res[[i]]$intercepts)
}

# set up a loop to simulate data and compute mean error
error_distributions <- list()  # store the error distributions here

## Simulate 100 models based on individual parameters
# This code is almost identical to the code by Haslbeck et al. (2021) 
# https://github.com/jmbh/tvvar_paper/blob/master/Tutorials/tutorial_mgm.R
time_before_sim <- Sys.time()
set.seed(2021)
for(i in pps){
  n <- sum(tvvar_res[[i]]$call$data_lagged$included) # number of included time points 
  for(j in 1:nIter){
    l_data[[j]] <- mlVAR::simulateVAR(pars = pars_sim[[i]], 
                                      means = intercepts_sim[[i]], 
                                      Nt = n, 
                                      residuals = 1)  
    # now fit a time-varying model on the data
    l_model[[j]] <- mgm::tvmvar(data = l_data[[j]],
                                type = rep("g", 6),
                                level = rep(1, 6), 
                                lambdaSel = "EBIC",                
                                lambdaGam = 0,
                                consec = 1:n,                  
                                estpoints = seq(0, 1, length = 20), 
                                bandwidth = df_min_bw[df_min_bw$id == i,]$bw_final,
                                lags = 1,
                                scale = TRUE,
                                pbar = FALSE, 
                                signInfo = FALSE)
    # store the prediction error and compute mean over all variables
    pred_res_sim <- predict(object = l_model[[j]], 
                        data = l_data[[j]], 
                        tvMethod = "weighted")  
    l_meanerror[j] <- mean(pred_res_sim$errors$RMSE)
  }
  error_distributions[[i]] <- l_meanerror
}
time_after_sim <- Sys.time() - time_before_sim

# Or directly load results
# load(file = paste0(data_dir, "/Processed/sim_error_distributions.RData"))

# Evaluation --------------------------------------------------------------
# Goals: 
# 1. Compare empirical error to %5 quantile of sampling distribution
# Example way to plot it: 
hist(error_distributions[[1]], xlim=c(0.9,1.1)) # sampling distribution under H0
abline(v = error_emp[[1]], col="red") # empirical error
quantile(error_emp[[1]], probs = 0.05) # get 0.05 quantile

## Write a loop for the evaluation
# Compare 5%-Quantile of sampling distribution to empirical error
# Create a data frame to store the results in
eval_df <- data.frame(
  id = rep(NA, length(pps)),
  empirical = rep(NA, length(pps)),
  sampling_quant = rep(NA, length(pps)),
  significant = rep(NA, length(pps)),
  min_sampling_quant = rep(NA, length(pps)),
  diff_mean = rep(NA, length(pps))
)
# Store: id, empirical error, 5% quantile, significant or not?
# And: difference to mean of sampling distribution

for(i in pps){
  eval_df[i,"id"] <- i     # store id
  eval_df[i, "empirical"] <- error_emp[[i]]   # store empirical error
  # get 5% quantile of sampling distribution under H0
  eval_df[i, "sampling_quant"] <- quantile(error_distributions[[i]], probs = 0.05) 
  eval_df[i, "significant"] <- if(eval_df[i, "empirical"] < eval_df[i, "sampling_quant"]){
    "yes"} else "no"
  # get minimum value of sampling distribution under H0 (just out of curiosity)
  eval_df[i, "min_sampling_quant"] <- min(error_distributions[[i]])
  # store difference between mean of sampling distribution and empirical error
  eval_df[i, "diff_mean"] <- mean(error_distributions[[i]])-error_emp[[i]]  
}

### Create a table with results
# As in Script 2., use the "niceTable"-function from GitHub
# https://raw.githubusercontent.com/RemPsyc/niceplots/master/niceTableFunction.R
# Instead of parsing it directly, it is adapted here to show 4 instead of 2 digits
source(here::here("Scripts", "Auxiliary Functions.R"))

# Create a table with results of the test and save it
# t_stationarity_test<- niceTable_moredigits(eval_df)
# save_as_docx(t_stationarity_test,
#              path = paste0(data_dir, "/Processed/Stationarity Tests.docx"))


# -------------------------------------------------------------------------
# Animate Plots -----------------------------------------------------------
# -------------------------------------------------------------------------
# Use package "animation" to create a GIF for networks of each participant
library(animation)
for(i in pps){
  animation::saveGIF({
    for(ep in c(1:20)){                                       # choose 3 timepoints
      qgraph(t(tvvar_res[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
             layout = "circle",                                    # choose layout of network
             labels = var_labels_n,                                # proper labelling
             edge.color = t(tvvar_res[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
             mar = rep(5, 4),                                      # margins of the graph
             vsize=11, esize=6, asize=5,                           # adjust variable size
             maximum = .5,                                         # max. effect for plotting relations
             pie = pred_res[[i]]$tverrors[[ep]][, 3],              # add prediction error to each node
             title = paste0("Estimation point = ", ep, ", ID = ", i),
             title.cex = 0.8)}
  }, movie.name = paste0("network_gif_", i, ".gif"), 
  ani.height = 800, ani.width = 800, interval = 1)                 # interval := time per image 
}

# Animate ID 2 seperately for a PPT presentation
for(i in 2){
  animation::saveGIF({
    for(ep in c(1:20)){                                       # choose 3 timepoints
      qgraph(t(tvvar_res[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
             layout = "circle",                                    # choose layout of network
             labels = var_labels_n,                                  # proper labelling
             edge.color = t(tvvar_res[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
             mar = rep(5, 4),                                      # margins of the graph
             vsize=15, esize=7, asize=6,                         # adjust variable size
             maximum = .5,                                         # max. effect for plotting relations
             pie = pred_res[[i]]$tverrors[[ep]][, 3],              # add prediction error to each node
             title = ep,
             title.cex = 3)}
  }, movie.name = "ppt_animation_id2.gif", interval = 0.5,
  ani.width = 800, ani.height = 800)
}


# -------------------------------------------------------------------------
# Session Info ------------------------------------------------------------
# -------------------------------------------------------------------------

# Save sessionInfo
# info <- sessionInfo()
# saveRDS(info, file = paste0(data_dir, "/Processed/SessionInfo.Rdata"))




