#-----------------------------------------------------------------------#
#                               Analysis script
# Temporal Dynamics of depressive Symptomatology: 
# An idiographic time series analysis
# 
# Created by: Björn Siepe
#                                   Part 1
#                               Preprocessing
#-----------------------------------------------------------------------#


# Description -----------------------------------------------------------
# In this script, the following steps are performed: 
# 1. Proper naming of variables
# 2. Cutting off irrelevant data points before & after data collection
# 3. Identifying breaks of data collection >7 days
# 4. Identifying longest consecutive data collection span without break
# 5. Removing study breaks & choosing the longest time series for each id
# 6. Searching for weeks with few data points at the end of a time series
# 7. Applying exclusion criteria
# 8. Descriptives of missing data
# 9. Missing data imputation with Kalman Filter



# -------------------------------------------------------------------------
# Loading data/packages & preparations ------------------------------------
# -------------------------------------------------------------------------
# Loading packages
library(readr)            # Reading in csv
library(dplyr)            # Data manipulation etc.
library(imputeTS)         # Kalman Imputation
library(naniar)           # Visualize/Summarize missing data
library(here)             # makes sourcing easier

# Create a data directory
# Other users need to adapt this to their own system
data_dir <- "Data"

# Create a figure directory
figure_dir <-"Output/"

# Read in file
daily_data <- readr::read_csv(paste0(data_dir,"/Raw/STEADY Tagesdaten.csv"))
str(daily_data)  # looks good
# Of note: 2 participants are not included in this file as they dropped out
# very early

# Only select columns that are relevant for data analysis
daily_data <- daily_data %>% 
  select(id, date, study_phase, YEAR, MONTH, DAY, NS_TST, PHQ2_1, PHQ2_2,
         rumination, socialize, socialize_val)

# Rename the columns to have concise styling
daily_data <- daily_data %>% 
  rename_with(tolower) %>%   # all lowercase
  rename(sleep = ns_tst,     # rename individually
         phq1 = phq2_1,
         phq2 = phq2_2,
         rumin = rumination,
         soc_quant = socialize,
         soc_qual = socialize_val)

## Change id to consecutive numbers to make looping easier
# by converting id to a factor with one unique level per id
# As a safety net: explicitly save numbers that are assigned to each individual
id_coding <- daily_data %>% 
  mutate(id_n = as.numeric(factor(id,levels=unique(id)))) %>% 
  select(id, id_n) %>% 
  rename(SDD_Code = id) %>%     # rename for later joining purposes
  distinct()

# Save id coding for later use
 # save(id_coding,
 #     file = paste0(data_dir,"/Processed/id_codes.RData"))
# write.csv2(id_coding,
#            file = paste0(data_dir, "/Processed/id_codes.csv"))


# Remove old id and proceed with new numeric one
daily_data <- daily_data %>% 
  mutate(id = as.numeric(factor(id,levels=unique(id)))) 


# -------------------------------------------------------------------------
# Trim data ---------------------------------------------------------------
# -------------------------------------------------------------------------
# Issue: there are many irrelevant days for each individual
# First, there are empty days before and after study collection
# Second, there are gaps between study collection periods

# First tackling the issue of empty days before and after
# To do so, find the first and the last non-NA observation per id and then
# trim the observations around it

## Try to get the row no. of first and last row that is not completely missing
# Important to first get the original row number stored, so that this is not
# affected by any filtering
# Save it as na_ind for later inspection
# Get row number for each id, then include all rows with at least 1 non-NA in 
# selected columns
na_ind <- daily_data %>% 
  group_by(id) %>% 
  mutate(row = row_number()) %>%   
  filter(if_any(c(phq1,phq2,sleep,rumin,soc_qual,soc_quant), 
                ~!is.na(.))) %>%    
  group_by(id) %>% 
  summarize(first_obs = min(row),
            last_obs = max(row))     # indicator of first and last row per id

# Now we cut off all rows before na_ind$first_obs and after na_ind$last_obs
# then safe this as daily cut
daily_cut <- daily_data %>% 
  group_by(id) %>% 
  mutate(row = row_number()) %>%     # get row number for each id 
  left_join(na_ind, by = "id") %>%   # attach the NA row indicators
  group_by(id) %>% 
  filter(row >= first_obs, row <= last_obs) %>%  # cut off irrelevant rows
  select(-c(first_obs, last_obs)) %>%     # delete irrelevant columns
  ungroup()

# Identify study breaks ---------------------------------------------------
## Use logical indicator to show if a row has any non-NA in the 6 relevant variables
# need to use "!" so that this turns out "TRUE" when there is any non-missing
daily_cut$valid <- !apply(daily_cut[,7:12], 1,
                          function(x) all(is.na(x)))  

# for easier operations, recode logical to integer (TRUE = 1)
daily_cut$valid <- as.integer(daily_cut$valid)


# I noticed that there are some periods for which either morning or evening data 
# are missing for a longer time span, but not all items are missing 
# (so valid is still = 1)
# -> also need to exclude these gaps by creating a new indicator for each variable
# Use a function that searches for consec. number of missings per item
# called "search_na_gap"
source(here::here("Scripts", "Auxiliary Functions.R"))

# # Apply this function to each item
# sleep_nas <- search_na_gap(daily_cut$sleep)
# phq1_nas <- search_na_gap(daily_cut$phq1)
# phq2_nas <- search_na_gap(daily_cut$phq2)
# rumin_nas <- search_na_gap(daily_cut$rumin)
# soc_quant_nas <- search_na_gap(daily_cut$soc_quant)
# soc_qual_nas <- search_na_gap(daily_cut$soc_qual)
# 
# # Bind these together
# item_nas <- as.data.frame(cbind(sleep_nas, phq1_nas, phq2_nas, 
#                                 rumin_nas, soc_quant_nas, soc_qual_nas))
# 
# # Add a column that indicates if one of the items gives out TRUE for each row
# # these rows would need to be excluded
# # ind = 1 if all items can be used
# item_na_ind <- item_nas %>% 
#   mutate(ind = ifelse(rowSums(item_nas) == 0, 1, 0))
# 
# # To reiterate the difference between valid and ind:
# # Example 1: valid == 1, ind == 0 -> at least one item is non-missing, but at least one item misses for more than 7 days
# 
# # attach this to the data frame
# daily_cut$ind <- item_na_ind$ind
# 
# # Include all rows for which ind == 1
# daily_break_prep <- daily_cut %>% 
#   filter(ind == 1)
# 


# New Version 30.04. ------------------------------------------------------
# Do the same thing as before, but group by user first
# otherwise, we might erroneously exclude rows
test <- daily_cut %>% 
  group_by(id) %>% 
  mutate(across(c(phq1,phq2,sleep,rumin,soc_qual,soc_quant),
                ~search_na_gap(.),
                .names = "{.col}_na")) %>% 
  ungroup() %>% 
  mutate(ind = ifelse(rowSums(select(.,ends_with("na"))) == 0, 1, 0))

test_break_prep <- test %>% 
  filter(ind == 1)

anti_join(test_break_prep, daily_break_prep)

# This adds exactly one day (the first day) to id 22
daily_break_prep <- test_break_prep

# Find longest time series per individual ---------------------------------
# Search consecutive original row numbers, as these can include short 
# spans of missing values that will be imputed later on
# To do so, subtract the previous row number from the current one
# If rows are consecutive, consec == 1 (1 day difference)
# If rows are non-consecutive, consec is equal to difference in days
daily_break_prep <- daily_break_prep %>% 
  group_by(id) %>% 
  mutate(consec = row-lag(row)) %>% 
  ungroup()

# 23 NAs in consec represent the first entry of each individual, so we
# can set this to 1 (as every first day is non-missing)
daily_break_prep$consec[is.na(daily_break_prep$consec)] <- 1

## Search for the longest streak of consec == 1 and use this time span
# for each individual
# Create a counter (consec_counter) that indicates the number of
# subsequent row numbers. A value of 100 e.g. means that this is the 100th
# subsequent observation, irrespective of missingness
# consec_grp gives an identifier to each run of consecutive days
daily_break_prep <- daily_break_prep %>% 
  group_by(id) %>% 
  mutate(consec_grp = with(rle(consec), rep(seq_along(lengths), lengths))) %>%
  group_by(id, consec_grp) %>% 
  mutate(consec_counter = seq_along(consec_grp)) %>%
  ungroup()

## Use consec_counter and consec_grp to identify the relevant
# time span (i.e. the longest) from our data and exclude the rest
# Compute the length of each consec_grp, i.e. the length of a run of 
# consecutive observations
# Then extract the group indicator for each id that has the longest run
consec_ind <- daily_break_prep %>% 
  group_by(id, consec_grp) %>% 
  summarize(max_consec = max(consec_counter)) %>% 
  slice(which.max(max_consec))

## With this indicator, cut the data again
# What remains is the study period with which we will be working with
# Filter for those consec_grp that are equal to consec_group.y,
# which is the maximum consecutive run per id found in consec_ind 
daily_break_cut <- daily_break_prep %>% 
  left_join(consec_ind, by = c("id")) %>% 
  group_by(id) %>% 
  filter(consec_grp.x == consec_grp.y) %>% 
  ungroup()

# Double check to see if it worked: the maximum of consec_counter
# must be equal to max_consec
# It worked!
daily_break_cut %>% 
  group_by(id) %>% 
  mutate(test_consec = max(consec_counter)) %>% 
  distinct(test_consec, max_consec)


# Exclude last weeks with bad data ----------------------------------------
# Goal: Exclude "bad weeks", i.e. 7-day intervals with less than 3 valid 
# data entries at the end of the time series
# Note: Used 7-day-intervals instead of weeks here
# Why at the end of the time series? We are willing to impute data collection
# breaks during the study, but do not want to use irregularly spaced data at
# end of data collection. For more specific reasoning, see Preregistration

## Cut the data set into bins of 7 days for each individual
# call the respective variable week_bin 
# This command repeats increasing integers in bins of 7 days
# seq_len(n()/2)) simply specifies the maximum number of week bins we calculate
# this just needs to be any number larger than n()/7
daily_break_cut <- daily_break_cut %>% 
  group_by(id) %>% 
  mutate(week_bin = as.factor(rep(seq_len(n()/2),  
                                  each = 7, length.out = n()))) %>% 
  ungroup()

## Calculate the number of rows that have at least 1 non-NA for each bin
# num_val = 5 means that there were 5 valid days in that 7 day period
# Also calculate the relative amount of valid observations per week ("rel_val")
# This is esp. important as some week_bins are shorter than 7 days
week_ind <- daily_break_cut %>% 
  group_by(id, week_bin) %>%
  summarize(rel_val = sum(valid)/length(week_bin),
            num_val = sum(valid)) %>% 
  ungroup()

## Find the last weeks per id (which is equal to the maximum bin number)
# slice_max returns the 5 highest values of week_bin, so the last 5 weeks
# Find weeks with less than 3 out of 7 valid observations 
# at the end of a time series
week_ind %>% 
  group_by(id) %>% 
  slice_max(week_bin, n = 5) %>% 
  filter(rel_val < 3/7)

# This does only give us week 43 of id 7, which is not at the end of the time
# time series:
week_ind %>% 
  filter(id == 7) %>% 
  slice_max(week_bin)

# -> No need to exclude any week at the end of data collection!

# -------------------------------------------------------------------------
# Missing Data ------------------------------------------------------------
# -------------------------------------------------------------------------
# Goals: Exclude participants with n<130 or >30% item-wise missings 
# Exclude bad sleep data
# Get Description of Missing/Present Data
# Then apply Kalman Filter to impute missing data

# Summarize percentage of NAs per ID per variable 
desc_miss <- daily_break_cut %>% 
  group_by(id) %>% 
  mutate(obs = n()) %>% 
  summarise(across(c(phq1,phq2,sleep,rumin,soc_qual,soc_quant),
                   ~sum(is.na(.)/obs),
                   .names = "{.col}_miss"))

## 1. Exclusion criterion: Item-wise missings
# Find maximum missing percentage 
# -> no one has more that 30% item-wise missings
apply(desc_miss, 2, max)

## 2. Exclusion criterion: >=130 valid observations
# Only include valid days
daily_break_cut %>% 
  filter(valid == 1) %>% 
  group_by(id) %>% 
  count()

# ID 23, ID 12 and ID 14 have way too few observations -> get excluded
# We then have a clean dataset, ready for imputation and descriptive stats
daily_clean <- daily_break_cut %>% 
  filter(!(id == 23)) %>% 
  filter(!(id == 14)) %>% 
  filter(!(id == 12))

## 3. Cell-specific criterion: Irregular sleep data
# Aim to exclude sleep <0 and sleep > 1440 (24h)
# There is no such data!
daily_clean %>% 
  filter(sleep <0 | sleep > 1440)


# Now compute non missing observations per ID and save
nonna_obs_per_id <- daily_clean %>% 
  filter(valid == 1) %>% 
  group_by(id) %>% 
  count() %>% 
  rename(days = n)

# Compute full number of observations per ID and save 
obs_per_id <- daily_clean %>% 
  group_by(id) %>% 
  count() %>% 
  rename(days = n)

# # Write to .csv
# write.csv2(obs_per_id,
#            file = paste0(data_dir,"/Processed/Observations per ID.csv"),
#            row.names = FALSE)
# #
# # Write to .Rdata
# save(obs_per_id,
#      file = paste0(data_dir, "/Processed/obs_per_id.Rdata"))

## Out of interest: check for instances of soc_quant = 0 and soc_qual !=0
# There are some examples of id = 2 who rated Social Quantity as 0
# and quality as non-zero, sometimes even as 100
daily_clean %>% 
  filter(soc_quant == 0 & soc_qual != 0)

# -------------------------------------------------------------------------
# Missing Data Descriptives -----------------------------------------------
# -------------------------------------------------------------------------
# Get NA summaries per ID with naniar package
na_summ <- daily_clean %>%
  group_by(id) %>%
  naniar::miss_var_summary() %>% 
  ungroup()

# Reformat to wide for descriptive table 
# obtain missing percentages per ID
na_pct <- na_summ %>%
  tidyr::pivot_wider(names_from = variable,
                     values_from = c(n_miss, pct_miss)) %>% 
  select(id, pct_miss_sleep, pct_miss_phq1, pct_miss_phq2, pct_miss_rumin,
         pct_miss_soc_quant, pct_miss_soc_qual) %>% 
  mutate(mean_miss = rowMeans(.[,2:7])) %>%   # get rowMeans for miss_pct across items
  mutate(across(-c(id),
                ~round(., digits = 2)))         # round to 2 digits for percentages

# Obtain overall average missingness percentage over all IDs and Items
na_summ %>%
  tidyr::pivot_wider(names_from = variable,
                     values_from = c(n_miss, pct_miss)) %>% 
  select(id, pct_miss_sleep, pct_miss_phq1, pct_miss_phq2, pct_miss_rumin,
         pct_miss_soc_quant, pct_miss_soc_qual) %>% 
  mutate(mean_miss = rowMeans(.[,2:7])) %>% 
  summarize(mean_miss_ave = mean(mean_miss))

# # Write to Rdata for later use
# save(na_pct,
#     file = paste0(data_dir, "/Processed/na_percentages.Rdata"))

# -------------------------------------------------------------------------
# Kalman Filter -----------------------------------------------------------
# -------------------------------------------------------------------------
# Goal: Apply Kalman Filter to impute all missing data
# Keep both imputed and non-imputed time series (for possible comparisons)
# Specify a max gap of 7 days in the na_kalman function
# this implies that 7 days of consecutive missings or more are not imputed
# (which should not be present anyway, just making sure)
set.seed(2021)
daily_imp <- daily_clean %>% 
  group_by(id) %>% 
  mutate(across(c(phq1,phq2, sleep,rumin, soc_qual, soc_quant),  # columns to be imputed
                ~imputeTS::na_kalman(.x, maxgap = 7),            # specify imputation 
                .names = "{.col}_imp")) %>%                      # add "_imp" to new column names
  ungroup()

# Get rid of the non-imputed and irrelevant columns
# rename the imputed columns back to original names
# Save this as the dataset that will be used for all analyses
data <- daily_imp %>% 
  select(id, date, year, month, day, sleep_imp, phq1_imp, phq2_imp,
         rumin_imp, soc_quant_imp, soc_qual_imp) %>% 
  rename(sleep = sleep_imp,
         phq1 = phq1_imp,
         phq2 = phq2_imp,
         rumin = rumin_imp, 
         soc_quant = soc_quant_imp,
         soc_qual = soc_qual_imp)

# Save these data for descriptives + analysis
# save(data,
#      file = paste0(data_dir, "/Processed/data.RData"))



