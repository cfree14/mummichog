

# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(TMB)
library(TMBhelper)
library(tidyverse)
library(ggsidekick)
library(lubridate)

# Directories
codedir <- "original_model_cleaned/code/helper_functions"
datadir <- "original_model_cleaned/data"
plotdir <- "original_model_cleaned/figures"

# Source helper functions
source(file.path(codedir, "for_TMB.R"))


# Build data
################################################################################

# Loop through datafiles
i <- 1
for(i in 1:4){
  
  # Read data
  data_orig <- read.data(file.path(datadir, paste0("crk", i, ".dat")))
  
  # Build data
  data <- list(ncreeks=data_orig$ncreeks %>% as.numeric(),
               nrel=data_orig$nrel %>% as.numeric(),
               nsites=data_orig$nsites %>% as.numeric(),
               nperiods=data_orig$nperiods %>% as.numeric(),
               ntraps=data_orig$ntraps %>% as.numeric(),
               count.mat=data_orig$count.mat, 
               distances=data_orig$times %>% as.numeric(),
               times=data_orig$times %>% as.numeric())
  
  # Export data
  saveRDS(data, file=file.path(datadir, paste0("crk", i, ".Rds")) )
  
}
