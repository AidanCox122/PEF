## INTERANNUAL TRENDS
# This code chunk takes the annual average of both marine bird and mammal
# abundance and environmental conditions within each of the six zones.
# This is meant to compare how persistent interannual trends are influenced
# by persistent differences in environmental conditions and aims to match
# the temporal resolution of the training data with the patterns we aim to assess:


# setup -------------------------------------------------------------------

library(tidyverse) 

env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone)

# add environmental data to each mbm obs. 
interannual_mbm_grid <- 
  env_grid %>% 
  # remove grid cells that do not fall witin zones (cannot be used to train)
  filter(!is.na(zone)) %>% 
  # select variables for model
  dplyr::select(Date, zone, bathy:salt) %>% 
  # find average value in each zone on each cruise date
  group_by(Date, zone) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>%  # n = 1644
  ungroup() %>% 
  right_join(mbm_data, by = c('Date', 'zone')) %>%
  # lots of NAs from years prior to 2017
  filter(!is.na(bathy)) %>% 
  mutate(year = lubridate::year(Date)) %>% 
  # calculate average conditions in each zone in each year
  group_by(year, zone, Species_code) %>% 
  summarize_if(is.numeric, mean, na.rm = T)


