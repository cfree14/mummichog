

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
codedir <- "new_model/code/helper_functions"
datadir <- "new_model/data"
plotdir <- "new_model/figures"

# Read data
data_tags_orig <- readxl::read_excel(file.path(datadir, "Data_Mummichog_tagging_2022.xlsx"), sheet="Initial_tagging", na="NA")
data_recaps_orig <- readxl::read_excel(file.path(datadir, "Data_Mummichog_tagging_2022.xlsx"), sheet="Recaptures", na="NA")


# Format data
################################################################################

# Format tag data
data_tags <- data_tags_orig %>% 
  # Rename
  janitor::clean_names() %>% 
  # Format date
  mutate(date=lubridate::ymd(date)) %>% 
  # Format site
  mutate(site=recode(site, 
                     "Raritan_Downstream"="Raritan-Downstream", 
                     "Raritan_Upstream"="Raritan-Upstream", 
                     "Passaic_Downstream"="Passaic-Downstream"))

# Inspect
str(data_tags)
table(data_tags$site)

# Format recap data
data_recaps <- data_recaps_orig %>% 
  # Rename
  janitor::clean_names() %>% 
  # Format date
  mutate(date=lubridate::ymd(date)) %>% 
  # Format site
  mutate(site=recode(site,
                     "Raritan_Downstream"="Raritan-Downstream",
                     "Raritan_Upstream"="Raritan-Upstream",
                     "Passaic"="Passaic-Downstream")) %>% 
  # Absolute value of distance
  mutate(distance_m=abs(distance_m))

# Inspect
str(data_recaps)
table(data_recaps$site)



# Format data
################################################################################

# Read data
# Datafile is a list with following:
# ncreeks = # of creeks
# nrel = # of released fish
# nsites = # of sites
# nperiods = # of time periods
# ntraps = # of traos
# count.mat = 
# distances = distances of sites from tagging location (vector length of nsites)
# times = days since time of tagging (vector length of nperiods)

# Sites
sites <- c("Raritan-Downstream", "Raritan-Upstream", "Passaic-Downstream")

# Loop through sites and build data
i <- 2
for(i in 1:length(sites)){
  
  # Site
  site_do <- sites[i]
  
  # Tagging data
  data_tags_do <- data_tags %>% 
    filter(site==site_do)
  tagging_date <- data_tags_do$date %>% unique()
    
  # Recaps data
  data_recaps_do <- data_recaps %>% 
    filter(site==site_do)
  
  # Site key
  site_key <- data_recaps_do %>% 
    dplyr::select(river, site, lat, lon, distance_m) %>% 
    unique() %>% 
    arrange(distance_m)
  
  # Period key
  period_key <- data_recaps_do %>% 
    dplyr::select(date) %>% 
    unique() %>% 
    arrange(date) %>% 
    mutate(days_post_release=as.numeric(date - tagging_date))
  
  # Format counts
  counts_df <- data_recaps_do %>% 
    # Simplify
    dplyr::select(date, distance_m, replicate, length_mm_tagged) %>% 
    # Summarize
    group_by(date, distance_m, replicate) %>% 
    summarize(ntagged=sum(!is.na(length_mm_tagged))) %>% 
    ungroup() %>% 
    # Spread
    mutate(replicate=paste0("trap", replicate)) %>% 
    spread(key="replicate", value="ntagged") %>% 
    # Arrange
    arrange(date, distance_m)
  
  # Convert to matric
  count_mat <- counts_df %>% 
    dplyr::select(-c(date, distance_m)) %>% 
    sapply(., as.numeric) %>% 
    unname()
  
  # Parameters
  nreleased <- nrow(data_tags_do) %>% as.numeric()
  nsites <- nrow(site_key) %>% as.numeric()
  distances_m <- site_key$distance_m %>% as.numeric()
  nperiods <- nrow(period_key) %>% as.numeric()
  days_since_tagging <- period_key$days_post_release %>% as.numeric()
  ntraps <- ncol(count_mat) %>% as.numeric()

  # Merge data
  out <- list(ncreeks = 1, # of creeks
              nrel = nreleased, # of released fish
              nsites = nsites, # of sites
              nperiods = nperiods,# of time periods
              ntraps = ntraps, # of traps
              count.mat = count_mat, # count matrix
              distances = distances_m, # distances of sites from tagging location (vector length of nsites)
              times = days_since_tagging) # days since time of tagging (vector length of nperiods))
  
  # Export data
  filename <- paste0(tolower(site_do), ".Rds")
  saveRDS(out, file=file.path(datadir, filename))
  write.csv(counts_df, file=file.path(datadir, paste0(tolower(site_do), ".csv")), row.names=F)
  
}








