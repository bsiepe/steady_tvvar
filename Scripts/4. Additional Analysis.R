#-----------------------------------------------------------------------#
#                               Analysis script
# Temporal Dynamics of depressive Symptomatology: 
# An idiographic time series analysis
# 
# Created by: Björn Siepe
#                                   Part 4
#                             Additional Analysis
#-----------------------------------------------------------------------#
# Description -------------------------------------------------------------
# This script contains supplemental analyses and analyses that were not
# preregistered and take a smaller place in the paper.

# Packages & Data ---------------------------------------------------------
library(rrapply)    # working with nested lists
library(tidyverse)  # for everything
library(mgm)        # tvvar and var models
library(readr)      # read csv files
library(runner)     # calculate running statistics
library(magick)     # create join images for gif creation
library(gifski)     # create gifs
library(qgraph)     # network plots

# Load necessary data
# we need: data, bw_res, tvvar_res, pred_res, pps, bw_res_small, bw_res_large
# Create a data directory
data_dir <- "Data"
# Create a figure directory
figure_dir <-"Output/"
# Or only load raw data
load(paste0(data_dir, "/Processed/data.RData"))
# Set up the loops by creating a vector with unique IDs 
pps <- unique(data$id)   

# -------------------------------------------------------------------------
# One-Standard-Error-Rule -------------------------------------------------
# -------------------------------------------------------------------------
# Here, I try to implement the 1-SE-Rule by Hastie et al. 
# to choose the largest bandwidth that is within 1 SE of the 
# bandwidth with the lowest RMSE overall
# This should prevent overfitting and prevent the choice of models
# that vary too much over time

# Here is a formula for the computation of SEs in CV:
# stat.cmu.edu/~ryantibs/datamining/lectures/18-val1.pdf
# according to these slides, we have to calculate the SE 
# by first calculating the mean error for a fold and then computing the SD
# of the mean
# The same is done here in row 299: 
# https://github.com/jmbh/ARVAR/blob/master/aux_functions.R

# Load bw_res 
# load(file = paste0(data_dir, "/Processed/bw_results.Rdata"))

# Work with a deeply nested list of CV errors, convert it to dataframe
# "fullErrorFolds" contains all RMSE for all folds
nfolds <- 5 # number of folds in cross validation
df_rmse_ses_new <- list()
for(i in pps){
 df_rmse_ses_new[[i]] <- as.data.frame(rrapply::rrapply(bw_res[[i]]$fullErrorFolds, how = "melt") %>% 
    pivot_wider(names_from = "L3") %>%          # get names of variable 
    select(-c(error_mean, error_time)) %>%      # get rid of unneccessary columns (error_mean is already rounded seemingly)
    rename(bw = L1, fold = L2) %>% 
    unnest(errors) %>%                          # unpack nested errors
    group_by(bw, fold) %>%                      # group within a fold
    summarize(mean_error_fold = mean(errors)) %>%   # compute mean error per fold
    group_by(bw) %>%          
    mutate(se = sd(mean_error_fold)/sqrt(nfolds),
           mean_error = mean(mean_error_fold)) %>%  # compute standard error across folds for each bw
    distinct(bw, .keep_all = TRUE) %>%   # works like summarize but keeps the other columns
    select(-mean_error_fold) %>%         # delete unnecessary fold-wise means
    ungroup() %>%   
    mutate(bw = as.numeric(recode(bw,
           "1" = 0.01,  # replace bandwidth index number with actual value
           "2" = 0.12,
           "3" = 0.23,
           "4" = 0.34,
           "5" = 0.45,
           "6" = 0.56,
           "7" = 0.67,
           "8" = 0.78,
           "9" = 0.89,
           "10" = 1.00)))) 
}


# 1 SER for small iteration of bw -----------------------------------------
# Load small bw results
# load(file = paste0(data_dir, "/Processed/bw_results_small.RData"))

nfolds <- 5 # number of folds in cross validation
df_rmse_ses_small <- list()
for(i in 2){
  df_rmse_ses_small[[i]] <- as.data.frame(rrapply::rrapply(bw_res_small[[i]]$fullErrorFolds, how = "melt") %>% 
                                         pivot_wider(names_from = "L3") %>%          # get names of variable 
                                         select(-c(error_mean, error_time)) %>%      # get rid of unneccessary columns (error_mean is already rounded seemingly)
                                         rename(bw = L1, fold = L2) %>% 
                                         unnest(errors) %>%                          # unpack nested errors
                                         group_by(bw, fold) %>% 
                                         summarize(mean_error_fold = mean(errors)) %>%   # compute mean error per fold
                                         group_by(bw) %>%  
                                         mutate(se = sd(mean_error_fold)/sqrt(nfolds),
                                                mean_error = mean(mean_error_fold)) %>%  # compute standard error across folds for each bw
                                         distinct(bw, .keep_all = TRUE) %>%   # works like summarize but keeps the other columns
                                         select(-mean_error_fold) %>%         # delete unneccessary fold-wise means
                                         ungroup() %>%   
                                         mutate(bw = as.numeric(recode(bw,
                                                                       "1" = 0.001,  # replace bandwidth index number with actual value
                                                                       "2" = 0.002,
                                                                       "3" = 0.003,
                                                                       "4" = 0.004,
                                                                       "5" = 0.005,
                                                                       "6" = 0.005,
                                                                       "7" = 0.006,
                                                                       "8" = 0.007,
                                                                       "9" = 0.008,
                                                                       "10" = 0.009)))) 
}



# Plot 1SER ---------------------------------------------------------------
# We plot the 1 SER for id = 11
# First obtain MSE + SE for the optimal bandwidth for plotting
mse_se_11 <- df_rmse_ses_new[[11]] %>% 
  filter(bw == 0.12) %>% 
  mutate(mse_se = mean_error + se) %>% 
  pull(mse_se)

# Plot for the supplement
for(i in 11){
  tmp_bw_plot_11 <- ggplot(df_rmse_ses_new[[i]], aes(x = bw))+
          geom_point(aes(y = mean_error))+
          geom_line(aes(y = mean_error))+
          geom_errorbar(aes(ymin = mean_error-se, # errorbars with SE
                            ymax = mean_error+se), alpha = 0.4)+
          # include shaded region to show 1SER
          geom_ribbon(aes(ymin = min(mean_error), ymax = mse_se_11), 
                      fill = "lightblue", color = "lightblue", 
                      alpha = 0.5, linetype = 0)+
          labs(# title = paste0("ID", i),   # only for iterations over multiple ids
               x = "Bandwidth", y = "Prediction Error")+
          # Styling
          theme_light()+
          theme(legend.text = element_text(size = 8),
          panel.grid.major.x = element_blank(), # remove gridlines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),       # remove borders
          axis.line = element_line(colour = "black"))
  ggsave(tmp_bw_plot_11, filename = paste0(figure_dir,"bw_ser_id_", i, ".svg"),
         device = "svg", width = 16, height = 9, units = "cm")
}


# Caluclate best bw with 1 SER --------------------------------------------
# Instead of finding the bandwidth that minimizes the MSE
# now search for the largest bandwidth that is within the smallest MSE+1 SE
l_min_bw_ser <- list()
for(i in pps){
  l_min_bw_ser[[i]] <- df_rmse_ses_new[[i]] %>% 
    # calculate cutoff as mean error of bw that minimizes mean_error + the respective standard error
    # looks ugly, sorry
    mutate(cutoff = df_rmse_ses_new[[i]]$mean_error[which.min(df_rmse_ses_new[[i]]$mean_error)]+
             df_rmse_ses_new[[i]]$se[which.min(df_rmse_ses_new[[i]]$mean_error)]) %>% 
    # search for largest bandwidth that is below or equal to the cutoff
    filter(mean_error<cutoff) %>%
    # get old bandwidth (that minimzes mean_error)
    # and new bandwidth (the largest bw below the cutoff)
    summarize(bw_old = df_rmse_ses_new[[i]]$bw[which.min(df_rmse_ses_new[[i]]$mean_error)],
              bw_new = max(bw)) %>% 
    # get id column
    mutate(id = i)
}
# This is the table used for supplement 5
min_bw_ser <- do.call(rbind, l_min_bw_ser)

# Compare the bandwidths of id = 2 seperately, because here  
# a different bw was chosen in the 2nd iteration of bw selection
mse_se_02 <- df_rmse_ses_small[[2]] %>% 
  filter(bw == 0.009) %>%
  # get mse + standard error for bw = 0.009
  mutate(mse_se = mean_error+se ) %>%
  pull(mse_se)
# now compare this to mse of the initially selected bw of 0.1

# For id = 2, we would select bw of 0.1 
min(bw_res[[2]]$meanError) < mse_se_02 

# make niceTable out of the results
# eval(parse("https://raw.githubusercontent.com/rempsyc/niceplots/e8712039bd3eea82933a8371fcf6001c244988c4/niceTableFunction.R", 
#            encoding = 'UTF-8'))
# t_min_bw_ser <- niceTable(min_bw_ser)
# save_as_docx(t_min_bw_ser, 
#              path = paste0(data_dir, "/Processed/1SER_Results.docx"))
# rm(t_min_bw_ser)


# -------------------------------------------------------------------------
# Participant 2: More estimation points -----------------------------------
# -------------------------------------------------------------------------
# Idea: As the chosen bw for participant 2 is very small,
# the number of estimation points might actually make a noticeable difference
# Here, we double it from 20 to 40
tvvar_res_est40 <- list()
# Set seed again
set.seed(2021)
for(i in 2){
  data_ind <- subset(data, id == i)                               # create individual data for i
  data_ind <- data_ind[daily_vars]                                # only select self-report variables
  tvvar_res_est40[[i]] <- mgm::tvmvar(data = data_ind,                         
                                type = rep("g", 6),                    # all gaussian variables
                                level = rep(1, 6),                     # all continuous variables
                                lambdaSel = "EBIC",                    # use EBIC for selection of regularization parameter
                                lambdaGam = 0,                         # set EBIC hyperparameter
                                lags = 1,                              # choose VAR(1) model
                                consec = 1:nrow(data_ind),             # useless? indicate which timepoints are consecutive
                                scale = TRUE,                          # scale all variables (important bc of regularization)
                                estpoints = seq(0, 1, length = 40),    # estimation points
                                bandwidth = df_min_bw[df_min_bw$id == i,]$bw_final, # use bandwidth that was selected before for each id
                                pbar = TRUE)                           # progress bar to show if it is worth to get another coffee
  
}

# Plot networks
for(i in 2){
  tvvar_res_est40[[i]]$edgecolor[, , , ][tvvar_res_est40[[i]]$edgecolor[, , , ] == "darkgreen"] <- c("darkblue")
}



# Save as graphic
pdf(file = paste0(figure_dir, "network_plots_40estpoints_02.pdf"))
for(i in 2){
  par(mfrow=c(7,3)) 
  for(tp in c(1:20)){                                       # choose 3 timepoints
    print(qgraph(t(tvvar_res_est40[[i]]$wadj[, , 1, tp]),       # transpose the weights matrix for qgraph default 
                 layout = "circle",                                    # choose layout of network
                 labels = var_labels_n,                                  # proper labelling
                 edge.color = t(tvvar_res_est40[[i]]$edgecolor[, , 1, tp]),  # again, transpose for qgraph and get edgecolor out of model object
                 mar = rep(5, 4),                                      # margins of the graph
                 vsize=14, esize=15, asize=13,                         # adjust variable size
                 maximum = 1,                                         # max. effect for plotting relations 
                 # pie = pred_res[[i]]$tverrors[[tp]][, 3],              # do not add prediction error here
                 title = paste0("Estimation point = ", tp, ", ID = ", i),        
                 title.cex=0.9))
    
  }
  par(mfrow=c(7,3))
  for(tp in c(21:40)){                                       # choose 3 timepoints
    print(qgraph(t(tvvar_res_est40[[i]]$wadj[, , 1, tp]),       # transpose the weights matrix for qgraph default 
                 layout = "circle",                                    # choose layout of network
                 labels = var_labels_n,                                  # proper labelling
                 edge.color = t(tvvar_res_est40[[i]]$edgecolor[, , 1, tp]),  # again, transpose for qgraph and get edgecolor out of model object
                 mar = rep(5, 4),                                      # margins of the graph
                 vsize=14, esize=15, asize=13,                         # adjust variable size
                 maximum = 1,                                         # max. effect for plotting relations 
                 # pie = pred_res[[i]]$tverrors[[tp]][, 3],              # do not add prediction error here
                 title = paste0("Estimation point = ", tp, ", ID = ", i),        
                 title.cex=0.9))
    
  }
  
}
dev.off()

### Gif of 40 models
# Animate ID 2 
for(i in 2){
  animation::saveGIF({
    for(ep in c(1:40)){                                       # choose timepoints
      qgraph(t(tvvar_res_est40[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
             layout = "circle",                                    # choose layout of network
             labels = var_labels_n,                                  # proper labelling
             edge.color = t(tvvar_res_est40[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
             mar = rep(5, 4),                                      # margins of the graph
             vsize=15, esize=7, asize=6,                         # adjust variable size
             maximum = .5,                                         # max. effect for plotting relations
             # pie = pred_res[[i]]$tverrors[[ep]][, 3],              # let pred error out here
             title = ep,
             title.cex = 3)}
  }, movie.name = "animation_est40_id2_1.gif", interval = 0.1,
  ani.width = 800, ani.height = 800)
}


### Estimating all 539 models
tvvar_res_est539 <- list()
# Set seed again
set.seed(2021)
for(i in 2){
  data_ind <- subset(data, id == i)                               # create individual data for i
  data_ind <- data_ind[daily_vars]                                # only select self-report variables
  tvvar_res_est539[[i]] <- mgm::tvmvar(data = data_ind,                         
                                       type = rep("g", 6),                    # all gaussian variables
                                       level = rep(1, 6),                     # all continuous variables
                                       lambdaSel = "EBIC",                    # use EBIC for selection of regularization parameter
                                       lambdaGam = 0,                         # set EBIC hyperparameter
                                       lags = 1,                              # choose VAR(1) model
                                       consec = 1:nrow(data_ind),             # indicate which timepoints are consecutive
                                       scale = TRUE,                          # scale all variables (important bc of regularization)
                                       estpoints = seq(0, 1, length = 539),    # estimation points
                                       bandwidth = df_min_bw[df_min_bw$id == i,]$bw_final, # use bandwidth that was selected before for each id
                                       pbar = TRUE)                           # progress bar to show if it is worth to get another coffee
  
}
# save(tvvar_res_est539, file = paste0(data_dir, "/Processed/tvvar_res_est539_id2.Rdata"))


# Plot networks
for(i in 2){
  tvvar_res_est539[[i]]$edgecolor[, , , ][tvvar_res_est539[[i]]$edgecolor[, , , ] == "darkgreen"] <- c("darkblue")
}

### Gif of 539 models
# Compilation of this gif failed for 539 models in one
# Therefore, do this differently: 
# Print 539 individual images, then combine them into one GIF
for(i in 2){
    for(ep in c(1:539)){
      # need to add zeroes to file names for correct ordering 
      if(ep <10){
        png(filename = paste0(figure_dir,"539nets/00", ep, "tvvar_est539.png"),
            width = 800, height = 800)
        qgraph(t(tvvar_res_est539[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
               layout = "circle",                                    # choose layout of network
               labels = var_labels_n,                                  # proper labelling
               edge.color = t(tvvar_res_est539[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
               mar = rep(5, 4),                                      # margins of the graph
               vsize=12, esize=5, asize=4,                         # adjust variable size
               maximum = .5,                                         # max. effect for plotting relations
               # pie = pred_res[[i]]$tverrors[[ep]][, 3],              # let pred error out here
               title = ep,
               title.cex = 3)
        dev.off()
      }
      if(ep >=10 & ep<100){
        png(filename = paste0(figure_dir,"539nets/0", ep, "tvvar_est539.png"),
            width = 800, height = 800)
        qgraph(t(tvvar_res_est539[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
               layout = "circle",                                    # choose layout of network
               labels = var_labels_n,                                  # proper labelling
               edge.color = t(tvvar_res_est539[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
               mar = rep(5, 4),                                      # margins of the graph
               vsize=12, esize=5, asize=4,                         # adjust variable size
               maximum = .5,                                         # max. effect for plotting relations
               # pie = pred_res[[i]]$tverrors[[ep]][, 3],              # let pred error out here
               title = ep,
               title.cex = 3)
        dev.off()
      }
      if(ep >= 100){
        png(filename = paste0(figure_dir,"539nets/", ep, "tvvar_est539.png"),
            width = 800, height = 800)
        qgraph(t(tvvar_res_est539[[i]]$wadj[, , 1, ep]),       # transpose the weights matrix for qgraph default 
               layout = "circle",                                    # choose layout of network
               labels = var_labels_n,                                  # proper labelling
               edge.color = t(tvvar_res_est539[[i]]$edgecolor[, , 1, ep]),  # again, transpose for qgraph and get edgecolor out of model object
               mar = rep(5, 4),                                      # margins of the graph
               vsize=12, esize=5, asize=4,                         # adjust variable size
               maximum = .5,                                         # max. effect for plotting relations
               # pie = pred_res[[i]]$tverrors[[ep]][, 3],              # let pred error out here
               title = ep,
               title.cex = 3)
        dev.off()
      }
      }
}

### Create a gif from these images
# List all 539 image names
# see: https://www.nagraj.net/notes/gifs-in-r/
imgs <- list.files(paste0(figure_dir, "539nets/"), full.names = TRUE)

# Read images
img_list <- lapply(imgs, magick::image_read)

# Join images together
img_joined <- magick::image_join(img_list)

# Write to Gif
magick::image_write_gif(image = img_joined,
                    path = paste0(figure_dir, "539nets/nets539_id2.gif"))


# Bootstrap Violin Plots --------------------------------------------------
# Obtain autoregressive effect of Anhedonia for ID = 7 over all estimation
# points
# Reading in bootstrap server results and tvvar results
# rel_res <- readRDS(file = paste0(data_dir, "/Processed/BootstrappingResults2.rds"))
# load(file = paste0(data_dir, "/Processed/tvvar_res.RData"))

boot_dist <- list()
for(i in 1:20){
  boot_dist[[i]] <- data.frame(
    boot_par = rel_res[[7]]$bootParameters[1,1,1,i,],
    ep = i)
}
df_boot_dist <- do.call(rbind, boot_dist)

# now obtain actual point prediction
df_point_est <- data.frame(est = tvvar_res[[7]]$wadj[1,1,1,],
                           ep = 1:20)
# attach to bootstrap dataframe
df_boot_dist <- df_boot_dist %>% 
  left_join(df_point_est, by = "ep") 

# Create violin plot with boxplot
tmp_boot_plot <- ggplot(df_boot_dist,aes(x = factor(ep), y = boot_par))+
  geom_violin(width = 1)+
  geom_boxplot(width = 0.15)+                # small boxplot
  theme_classic()+
  # add point estimate
  geom_point(data = df_boot_dist, 
             aes(x = factor(ep), y = est), 
             col = "red", size = 1.5)+
  labs(x = "Estimation Point", y = "Parameter Value")
ggsave(tmp_boot_plot, filename = paste0(figure_dir,"boot_samp_violin.svg"),
       device = "svg", width = 16, height = 9, units = "cm")
rm(tmp_boot_plot)


