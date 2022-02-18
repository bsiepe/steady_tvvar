# README 

## Information
This file contains information about the data and other material that supplements the thesis "Temporal Dynamics of Depressive Symptomatology: An idiographic time series analysis".
The folder structure is as follows:

## Data
This folder contains all raw and processed data. Subdirectories include:
### Raw
This folder contains all raw data, more specifically: 
- **STEADY_Tagesdaten.csv**: Raw daily self-report data used for this study
- **Selbstbericht_Pre_Post.csv**: Pre-and post questionnaires, including the IDS-C used in the manuscript. 
- **Daten_STEADY_final.csv**:Demographic information and screening data.

### Processed
This folder contains all processed data as well as some R-Output. 
- **data.RData**: Final dataset with imputated data created by the script **1. Preprocessing.R**. 
- **id_codes.RData**: Id codes to match STEADY codes to my numerical codes.
- **bw_results.RData**: Rresults of first iteration of bandwidth selection. 
- **bw_results_small.RData**: Results of second iteration of bandwidth selection for small bandwidths.
- **bw_results_large.RData**: Results of second iteration of bandwidth selection for large bandwidths.
- **tvvar_res.RData**: Resuts of time-varying vectorautoregression.
- **sim_error_distributions.RData**: Results of simulations for hypothesis test. 
- **SessionInfo.RData**: Session Info 

#### Output
This folder contains various Figures and animations used in the manuscript, supplement or not used yet. 
Importantly, **nets539_id2.gif** is the GIF of all 539 models estimated for participant 2 that is referenced in a footnote of the manuscript.

## Script
This folder contains an R-Project as well as all scripts needed to reproduce the analyses in the study. The scripts are numbered and should be used in that sequence. Still, each script contains loading comments to load all relevant data or functions from other scripts. These are the scripts:
- **1. Preprocessing.R**: Contains missing data procedures and imputation, creates final data set. 
- **2. Descriptives.R**: Contains all descriptives used in the manuscript. 
- **3. Analysis.R**: Contains all main and preregistered analyses in the manuscript. 
- **4. Additional Analyses.R**: Contains all analyses for the supplementary material as well as further visualizations.
- **5. Block Bootstrap Server.R**: Contains code to perform the bootstrap stability analyses for computation on a Linux server. 
- **Auxiliary Functions.R**: Contains aux functions that are sourced throughout the other scripts.

## Manuscript
This folder contains the manuscript as well as the supplemental materials.
