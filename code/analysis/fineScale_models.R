
# setup -------------------------------------------------------------------

library(tidyverse) 
library(PerformanceAnalytics)
library(mgcv)

source('code/functions.R')


env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone)

