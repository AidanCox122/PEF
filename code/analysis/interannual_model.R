## INTERANNUAL TRENDS
# This code chunk takes the annual average of both marine bird and mammal
# abundance and environmental conditions within each of the six zones.
# This is meant to compare how persistent interannual trends are influenced
# by persistent differences in environmental conditions and aims to match
# the temporal resolution of the training data with the patterns we aim to assess:


# setup -------------------------------------------------------------------

library(tidyverse) 
library(PerformanceAnalytics)
library(mgcv)

source('code/functions.R')

env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone) %>% 
  mutate(zone = factor(zone))

# add environmental data to each mbm obs. 
interannual_mbm_grid <- 
  env_grid %>% 
  # remove grid cells that do not fall witin zones (cannot be used to train)
  filter(!is.na(zone)) %>% 
  # change zone to a factor
  mutate(zone = factor(zone)) %>% 
  # select variables for model
  dplyr::select(Date, zone, bathy:salt) %>% 
  # scale predictors
  mutate_if(is.numeric, base::scale) %>% 
  # find average value in each zone on each cruise date
  group_by(Date, zone) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>%  # n = 1644
  ungroup() %>% 
  right_join(mbm_data, by = c('Date', 'zone')) %>%
  # lots of NAs from years prior to 2017
  filter(!is.na(bathy)) %>% # n = 672
  mutate(year = lubridate::year(Date)) %>% 
  # calculate average conditions in each zone in each year
  group_by(year, zone, Species_code) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>% 
  ungroup() %>% 
  mutate(countInt = round(Count,0))

# scale the variables 

# model training -----------------------------------------------------------
 
interannual_mbm_grid %>% 
  dplyr::select('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt') %>% 
  chart.Correlation()

## GL ----------------------------------------------------------------------

# forward selection 1
getWeightIANN(test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'GL',
              training = interannual_mbm_grid)
# bathymetry is the best predictor

# forward selection 2
getWeightIANN(base = c('bathy'), test = c('topog', 'dist', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'GL',
              training = interannual_mbm_grid)
# distance from shore is the best predictor

# forward selection 3
getWeightIANN(base = c('bathy', 'dist'),
              test = c('phyto', 'sst', 'temp_sd', 'salt'),
              species = 'GL',
              training = interannual_mbm_grid)
# salinity is the best predictor

# forward selection 4
getWeightIANN(base = c('bathy', 'dist', 'salt'),
              test = c('phyto', 'sst', 'temp_sd'),
              species = 'GL',
              training = interannual_mbm_grid)
# phytoplankton is the best predictor

# forward selection 5
getWeightIANN(base = c('bathy', 'dist', 'salt', 'phyto'),
              test = c('sst', 'temp_sd'),
              species = 'GL',
              training = interannual_mbm_grid)
# temperature standard deviation is the best predictor

# forward selection 6
getWeightIANN(base = c('bathy', 'dist', 'salt', 'phyto', 'temp_sd'),
              test = c('sst'),
              species = 'GL',
              training = interannual_mbm_grid)

# roughly 50/50 w. base model so stopping here

### best model predicting interannual variability in GL: -----
GL_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=4)+s(dist,k=4)+s(salt,k=4)+s(phyto,k=4)+s(temp_sd,k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'GL')))
summary(GL_interann_mod)
# r-squ. adj = 0.619

## CoMu ----------------------------------------------------------------------
# forward selection 1
getWeightIANN(test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# distance from shore is the best predictor

# forward selection 2
getWeightIANN(base = c('dist'),
              test = c('bathy', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# phytoplankton concentration is the best predictor

# forward selection 3
getWeightIANN(base = c('dist', 'phyto'),
              test = c('bathy', 'sst', 'temp_sd', 'salt'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# salinity is the best predictor

# forward selection 4
getWeightIANN(base = c('dist', 'phyto', 'salt'),
              test = c('bathy', 'sst', 'temp_sd'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# sst is the best predictor

# forward selection 5
getWeightIANN(base = c('dist', 'phyto', 'salt', 'sst'),
              test = c('bathy', 'temp_sd'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# sst_sd is the best predictor

# forward selection 6
getWeightIANN(base = c('dist', 'phyto', 'salt', 'sst', 'temp_sd'),
              test = c('bathy'),
              species = 'CoMu',
              training = interannual_mbm_grid)
# bathymetry is the best predictor 
# BUT SST loses signif when bathy included so stoping at FS 5

### best model for CoMu ----
CoMu_interann_mod <- 
  gam(formula = countInt~s(dist,k=4)+s(phyto,k=4)+s(salt,k=4)+s(sst,k=4)+s(temp_sd,k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = data)
summary(CoMu_interann_mod)
# r-sq. adj. = 0.468

## HSeal ----------------------------------------------------------------------

# forward selection 1
getWeightIANN(test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'HSeal',
              training = interannual_mbm_grid)
# bathymetry is the best predictor

# forward selection 2
getWeightIANN(base = c('bathy'),
              test = c('topog', 'dist', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'HSeal',
              training = interannual_mbm_grid)
# phytoplankton concentration is the best predictor

# forward selection 3
getWeightIANN(base = c('bathy', 'phyto'),
              test = c('topog', 'dist', 'sst', 'temp_sd', 'salt'),
              species = 'HSeal',
              training = interannual_mbm_grid)
# topography is next best predictor

# forward selection 4
getWeightIANN(base = c('bathy', 'phyto', 'topog'),
              test = c('sst', 'temp_sd', 'salt'),
              species = 'HSeal',
              training = interannual_mbm_grid)
# sst is next best predictor

# forward selection 5
getWeightIANN(base = c('bathy', 'phyto', 'topog', 'sst'),
              test = c('temp_sd', 'salt'),
              species = 'HSeal',
              training = interannual_mbm_grid)
# null model is next best predictor

### best model for HSeals ----
HSeal_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=4)+s(phyto,k=4)+s(topog,k=4)+s(sst,k=4),
    family = poisson,
    offset = log(Effort_sqkm),
    data = data)
summary(HSeal_interann_mod)
# r-squ. adj. = 0.522


## HPorp ----------------------------------------------------------------------
# forward selection 1
getWeightIANN(test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'HPorp',
              training = interannual_mbm_grid)
# bathymetry is the best predictor

# forward selection 2
getWeightIANN(base = c('bathy'),
              test = c('topog', 'dist', 'phyto', 'sst', 'temp_sd', 'salt'),
              species = 'HPorp',
              training = interannual_mbm_grid)
# sst is the best predictor

# forward selection 3
getWeightIANN(base = c('bathy', 'sst'),
              test = c('topog', 'dist', 'phyto', 'temp_sd', 'salt'),
              species = 'HPorp',
              training = interannual_mbm_grid)
# salinity is the best predictor
# BUT salt makes sst insignifiant and so we'll stop here

### best model for HSPorps ----
HPorp_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=4) + s(sst, k=4) + s(salt, k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = data)
summary(HPorp_interann_mod) # this is significant!


# Leave-one-out Validation ------------------------------------------------


