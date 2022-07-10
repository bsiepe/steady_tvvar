#-----------------------------------------------------------------------#
#                               Analysis script
# Temporal Dynamics of depressive Symptomatology: 
# An idiographic time series analysis
# 
# Created by: Björn Siepe
#                                   Part 2
#                               Descriptives
#-----------------------------------------------------------------------#

# Description -------------------------------------------------------------
# In this script, the following steps are performed: 
# 1. Describe individual items
# 2. Compute IDS-C Pre and Post 
# 3. Get descriptives for sample
# 4. Build Histogram and Time Series plots for all vars

# -------------------------------------------------------------------------
# Loading data/packages & preparations ------------------------------------
# -------------------------------------------------------------------------
library(dplyr)              # data manipulation
library(moments)            # skewness and kurtosis
library(readxl)             # reading Excel files
library(tidyr)              # data reshaping
library(gridExtra)          # grid graphics arrangement
library(ggplot2)            # ...
library(here)               # makes sourcing easier

# Parse niceTable() function from Github for creating nice tables in APA style
eval(parse("https://raw.githubusercontent.com/rempsyc/niceplots/e8712039bd3eea82933a8371fcf6001c244988c4/niceTableFunction.R", 
           encoding = 'UTF-8'))

# Other users need to adapt this to their own system
data_dir <- "Data"

# Create a figure directory
figure_dir <-"Output/"

# Load processed data
load(paste0(data_dir, "/Processed/data.RData"))

# -------------------------------------------------------------------------
# Individual descriptives of daily items ----------------------------------
# -------------------------------------------------------------------------
# Compute person-specific means, SDs, skew and kurtosis for every participant
desc_summ <- data %>% 
  group_by(id) %>% 
  summarise(across(c(phq1,phq2, sleep,rumin, soc_quant, soc_qual),
                   list(mean = mean, sd = sd, 
                        skew = moments::skewness, kurt = moments::kurtosis), 
                   .names = "{.col}_{.fn}")) # proper column naming


## Compute frequency count of individual responses
# This can highlight floor/ceiling effects
# Use the auxiliary function "floorceiling" to compute relative frequencies
source(here::here("Scripts", "Auxiliary Functions.R"))

# For example, find relative frequencies for phq1
floorceiling(data, phq1)

# -------------------------------------------------------------------------
# Pre-Post Data -----------------------------------------------------------
# -------------------------------------------------------------------------
# Read in Pre-Post Questionnaires for using the IDS
pre_post <- readxl::read_excel(paste0(data_dir, "/Raw/Selbstbericht_Pre_Post.xlsx"), 
                               sheet = "Daten")

## Remove superficial columns
# only use information of the IDS
pre_post <- pre_post %>% 
  select(SDD_Code, Anmerkungen, Messpunkt, Datum,
         IDSC_I01:IDSC_I30) %>% 
  mutate(SDD_Code = as.factor(SDD_Code)) %>%      # for later manipulation
  mutate(Datum = as.Date(Datum))
## Recode the SDD_Code to comply with id coding in daily dataset
# load the id codes generated in the Preprocessing script
load(paste0(data_dir,"/Processed/id_codes.RData"))

# Recode the ids to include an underscore and additional 0
# to comply with coding in pre_post
id_coding$SDD_Code <- as.factor(gsub("SDY", "SDY_0", id_coding$SDD_Code))


## Perform recoding and exclude irrelevant ids
# these ids are not used in this study as they dropped out early
# and do not have daily data
# This removes SDY_002 and SDY_006, who did not participate
pre_post <- pre_post %>% 
  left_join(id_coding, by = "SDD_Code") %>% 
  mutate(id = id_n) %>%         # id_n contains correct numerical id
  select(-id_n) %>%             # remove irrelevant column
  select(id, SDD_Code, everything()) %>%   # reorder columns
  filter(!is.na(id)) %>%        # remove those that are not included in daily data
  filter(!is.na(Datum)) %>%     # remove irrelevant placholder rows + canceled days
  filter(if_any(c(IDSC_I01:IDSC_I30), # remove rows with only missings
              ~!is.na(.)))
## Recode Missings 
# As per the codebook:
# 999 = missing
# 998 = does not apply
# We treat both as missing
pre_post[pre_post[,] == 999] <- NA
pre_post[pre_post[,] == 998] <- NA
  
## Find data points close to daily data phase for each id
# Issue: As some daily data is excluded, the first IDS questionnaire might not be closest
# to the start of the daily data phase for everyone
# Solution: Search for nearest date before daily data starts and after it ends

## First: find out start and end dates of daily data for each id
# load(paste0(data_dir, "Processed/data.Rdata"))
first_last <- data %>% 
  group_by(id) %>% 
  summarize(first_date = min(date),
         last_date = max(date))

# Attach the dates to the pre_post data frame
pre_post <- pre_post %>% 
              left_join(first_last,
              by = "id")

# Find nearest IDS before and after daily data period
# Find pre-questionnaire by minimizing absolute distance to data collection
# Then build a sum score for IDSC
idsc_pre <- pre_post %>% 
  group_by(id) %>% 
  filter(Datum <= first_date) %>%        # select only dates before daily data started
  mutate(pre_date_diff = as.numeric(first_date-Datum)) %>% 
  slice(which.min(abs(pre_date_diff))) %>% # select date minimizing difference
  ungroup %>% 
  mutate(idsc_pre_sum = rowSums(select(.,starts_with("IDSC")), na.rm = TRUE))

# Find post-questionnaire
idsc_post <- pre_post %>% 
       group_by(id) %>% 
       filter(Datum >= first_date) %>%        # select only dates after daily data ended
       mutate(post_date_diff = as.numeric(last_date-Datum)) %>% 
       slice(which.min(abs(post_date_diff))) %>% 
       ungroup %>% 
       select(-c(IDSC_I09B, IDSC_I09C)) %>%   # Exclude both non-Likert items
       mutate(idsc_post_sum = rowSums(select(.,starts_with("IDSC")), na.rm = TRUE))

# Build new data.frame
idsc_sum <- idsc_pre %>% 
  left_join(idsc_post %>% select(id, idsc_post_sum), by = "id") %>% 
  select(id, idsc_pre_sum, idsc_post_sum)

# # save for later use
# save(idsc_sum,
#     file = paste0(data_dir, "/Processed/idsc_sumscores.Rdata"))

# -------------------------------------------------------------------------
# Sample Description Table ------------------------------------------------
# -------------------------------------------------------------------------
# Read in description of sample
screen <- readxl::read_excel(paste0(data_dir, "/Raw/Daten_STEADY_final.xlsx"), 
                               sheet = "Screenings")

## Load information about missingness and data collection from
# 1. Preprocessing
# Load Missing Info per ID
load(paste0(data_dir, "/Processed/na_percentages.Rdata"))

# Load number of observations per ID
load(paste0(data_dir, "/Processed/obs_per_id.Rdata"))

# Rename identification variable
colnames(screen)[1] <- "SDD_Code"

# Attach our IDs
screen <- screen %>% 
  left_join(id_coding, by = "SDD_Code") %>% 
  mutate(id = id_n) %>%         # id_n contains correct numerical id
  select(-id_n) %>%             # remove irrelevant column
  select(id, SDD_Code, everything()) %>%   # reorder columns for clarity
  filter(!is.na(id)) %>%        # remove those not included in daily data
  filter(!(id == 23)) %>%       # excluded participants (see 1. Preprocessing)
  filter(!(id == 14)) %>% 
  filter(!(id == 12))

## Age summary 
mean(screen$Age)
sd(screen$Age)
range(screen$Age)

## Sex summary  
# (Female = 1)
screen %>% 
  count(Sex)

## Nationality summary 
# 1 = German, 2 = dual citizen, 3 = non-German
screen %>% 
  count(Nat)

## Education summary
# 8 = have passed their Abitur
screen %>% 
  count(School_01)

## "Screen" is rich in information (e.g. detailed SCID)
# here, shorten it for summary table
# Variable names can be accessed in the document:
# "Variablen_Itemnamen_Codierungen.xlsx"
screen_sh <- screen %>% 
  select(id, Sex, Age) %>% 
  mutate(Sex = as.factor(Sex),
         Age = as.factor(Age)) %>%
  left_join(obs_per_id, by = "id") %>%             # attach observations per ID 
  left_join(na_pct %>% select(id, mean_miss), 
            by = "id") %>%                         # attach mean NA % per ID
  left_join(idsc_sum, by = "id") %>%               # attach IDSC sum scores
  mutate(Sex = recode(Sex, "1" = "F", "2" = "M")) %>%   # recode Sex
  mutate(id = as.factor(id))                       # mutate id after left_joining to avoid issues with different formats

# Work with this information in Word
summ_table <- niceTable(screen_sh)
save_as_docx(summ_table, 
             path = paste0(data_dir, "/Processed/Summary Table.docx"))

## Create a dataframe with treatment information
treat_sh <- screen %>% 
  select(id, SDD_Code, Med_01_02, Med_02_02, Med_03_02, Med_04_02, Med_05_02,
         PT_01_02, PT_02_02, PT_03_02, PT_04_02)

## Create a dataframe with employment information
# see file "Variablen_Itemnamen_Codierungen.xlsx"
employ_sh <- screen %>% 
  select(id, SDD_Code, Work_01_01:Work_02) %>% 
  rename(fulltime = Work_01_01,
         parttime = Work_01_02,
         selfempl = Work_01_03,
         partialretire = Work_01_04,
         minimalemploy = Work_01_05,
         eurojob = Work_01_06,
         irregular = Work_01_07,
         search = Work_01_08,
         training = Work_01_09,
         reeducation = Work_01_10,
         milcivservice = Work_01_11,
         parentleave = Work_01_12,
         retired = Work_01_13,
         school = Work_01_14,
         permunable = Work_01_15,
         housemanwoman = Work_01_16,
         other = Work_01_17,
         weeklyhours = Work_02)  # rename according to proper meaning

## Create a dataframe with diagnoses
diag_sh <- screen %>%
  select(id, SDD_Code, Erk_00:Erk_02_02)


### Create table for missingness for the supplement with na_pct
# reorder some variables
na_pct <- na_pct %>% 
  select(id, pct_miss_phq1, pct_miss_phq2, pct_miss_sleep, 
         pct_miss_rumin, pct_miss_soc_quant, pct_miss_soc_qual)

# t_supp_miss <- niceTable(na_pct)
# save_as_docx(t_supp_miss, 
             # path = paste0(data_dir, "/Processed/Missing Data Supplementary.docx"))

# -------------------------------------------------------------------------
# Item Histograms ---------------------------------------------------------
# -------------------------------------------------------------------------
# Goal: Create histograms for each variable for each ID
# Use grid.arrange to distribute plots over multiple pages
pdf(paste0(data_dir, "/Processed/histograms_allvars.pdf"))
gridExtra::grid.arrange(
  ggplot(data = data)+
  geom_histogram(aes(x = phq1), binwidth = 1)+
  facet_wrap(.~id, ncol = 5),
  ggplot(data = data)+
    geom_histogram(aes(x = phq2), binwidth = 1)+
    facet_wrap(.~id, ncol = 5), nrow = 2)
gridExtra::grid.arrange(ggplot(data = data)+
    geom_histogram(aes(x = sleep), binwidth = 20)+
    facet_wrap(.~id, ncol = 5),
  ggplot(data = data)+
    geom_histogram(aes(x = rumin), binwidth = 1)+
    facet_wrap(.~id, ncol = 5), nrow = 2)
# adjust binwidth to larger scale of social quant + qual
gridExtra::grid.arrange(ggplot(data = data)+
    geom_histogram(aes(x = soc_quant), binwidth = 10)+
    facet_wrap(.~id, ncol = 5),
  ggplot(data = data)+
    geom_histogram(aes(x = soc_qual), binwidth = 10)+
    facet_wrap(.~id, ncol = 5), nrow = 2)
dev.off()

# -------------------------------------------------------------------------
# Item Time Series --------------------------------------------------------
# -------------------------------------------------------------------------
# Create a time series plot for every variable to show development over time
# First, need to add a study_day variable that gives day per id
data_plot <- data %>% 
  group_by(id) %>% 
  mutate(study_day = row_number())

pdf(paste0(data_dir, "/Processed/timeseries_allvars.pdf"))
ggplot(data = data_plot)+
    geom_line(aes(x = study_day, y = phq1), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
    theme_classic()
ggplot(data = data_plot)+
    geom_line(aes(x = study_day, y = phq2), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
    theme_classic()
ggplot(data = data_plot)+
    geom_line(aes(x = study_day, y = sleep), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
    theme_classic()
ggplot(data = data_plot)+
    geom_line(aes(x = study_day,y = rumin), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
    theme_classic()
ggplot(data = data_plot)+
    geom_line(aes(x = study_day,y = soc_quant), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
  theme_classic()
ggplot(data = data_plot)+
    geom_line(aes(x = study_day,y = soc_qual), size = 0.1, colour = '#0d0887')+
    facet_wrap(.~id, ncol = 5)+
  theme_classic()
dev.off()

# remove data_plot for clearness
rm(data_plot)
