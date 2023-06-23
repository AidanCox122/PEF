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
      data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu')))
summary(CoMu_interann_mod)
# r-sq. adj. = 0.795

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
# topography is next best predictor but not statistically significant stop here

### best model for HSeals ----
HSeal_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=4)+s(phyto,k=4),
    family = poisson,
    offset = log(Effort_sqkm),
    data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal')))
summary(HSeal_interann_mod)
# r-squ. adj. = 0.617


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

# forward selection 4
getWeightIANN(base = c('bathy', 'sst', 'salt'),
              test = c('topog', 'dist', 'phyto', 'temp_sd'),
              species = 'HPorp',
              training = interannual_mbm_grid)
# null model is best

### best model for HSPorps ----
HPorp_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=4)+s(sst, k=4)+s(salt, k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp')))
summary(HPorp_interann_mod) # this is (nearly) significant!
# adj. r-squared: 0.621

# Leave-one-out Validation ------------------------------------------------
# This chunk performs leave-one-year-out validation on course-scale habitat models:

# pivot interannual model to wide format so MBM columns can be used for different models w.o filtering
cv_base <- 
  interannual_mbm_grid %>% 
  dplyr::select(year, zone, Effort_sqkm, Species_code, bathy:salt, countInt) %>% 
  pivot_wider(names_from = 'Species_code', values_from = 'countInt')

# loop that filters out each year, trains each species model, and computes/stores validation metrics
LYO_results <- data.frame()
LYO_raw <- data.frame()
for (y in unique(cv_base$year)) {
  ## create train and test sets
  counter = 1
  train <- cv_base %>% filter(year != y)
  test <- cv_base %>% filter(year == y)
  ## train a model for each species
  lm.HP <- gam(formula = HPorp~ s(bathy,k=4)+s(sst, k=4)+s(salt, k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.HS <- gam(formula = HSeal~ s(bathy,k=4)+s(phyto,k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.CM <- gam(formula = CoMu~s(dist,k=4)+s(phyto,k=4)+s(salt,k=4)+s(sst,k=4)+s(temp_sd,k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.GL <- gam(formula = GL~ s(bathy,k=4)+s(dist,k=4)+s(salt,k=4)+s(phyto,k=4)+s(temp_sd,k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  ## apply these models to test data
  test.full <- 
    test %>% 
    mutate(
      # predict values for harbor porpoise
      HP.p = predict(lm.HP, newdata = test, type = "response", se = FALSE),
      # round non-integer values produced be predict() function
      HP.p = round(HP.p, digits = 0),
      # repeat for Harbor seal
      HS.p = predict(lm.HS, newdata = test, type = "response", se = FALSE),
      HS.p = round(HS.p, digits = 0),
      # Common Murre
      CM.p = predict(lm.CM, newdata = test, type = "response", se = FALSE),
      CM.p = round(CM.p, digits = 0),
      # Glaucous-winged gull
      GL.p = predict(lm.GL, newdata = test, type = "response", se = FALSE),
      GL.p = round(GL.p, digits = 0)) %>% 
    # convert zero values to very small numbers to avoid infinite predition error
    mutate_all(across(everything(), ~ if_else(
      . == 0,
      0.000001,
      .)))
    # calculate prediction error for each species in each zone
    mutate(
      d.HP = (((HPorp - HP.p) / HPorp) * 100),
      d.HS = (((HSeal - HS.p) / HSeal) * 100),
      d.CM = (((CoMu - CM.p) / CoMu) * 100),
      d.GL = (((GL - GL.p) / GL) * 100)
    )
  
  ## compute and store validation metrics
  # prediction error and adj. r-squared
  test$HS[test$HS == 0] <- 0.001 # prevent infinite prediction error
  test$HP[test$HP == 0] <- 0.001 # prevent infinite prediction error
  test$CM[test$CM == 0] <- 0.001 # prevent infinite prediction error
  test$GL[test$GL == 0] <- 0.001 # prevent infinite prediction error
  test <- test %>% mutate(
    
  )
  data <- data.frame(
    year = y,
    species = rc("HP", "HS", "CM", "GL"),
    R.squared = c(summary(lm.HP)$r.sq, summary(lm.HS)$r.sq, summary(lm.CM)$r.sq, summary(lm.GL)$r.sq),
    Dev.Expl = c(summary(lm.HP)$dev.ex, summary(lm.HS)$dev.ex, summary(lm.CM)$dev.ex, summary(lm.GL)$dev.ex)
  )
  
  LYO_raw <- rbind(LYO_raw, test[,c(1,2,13:20)])
  LYO_results <- rbind(LYO_results, data)
  print(y)
  rm(data, train, test, lm.HP, lm.HS, lm.CM, lm.GL)
}

# same process but for zones
LZO_results <- data.frame()
LZO_raw <- data.frame()
for (z in unique(cv_base$Zone)) {
  ## create train and test sets
  train <- cv_base %>% filter(Zone != z)
  test <- cv_base %>% filter(Zone == z)
  ## train a model for each species
  lm.HP <- gam(HP ~ s(phyto, k=4) , data = train)
  lm.HS <- gam(HS ~ s(bathy,k=3) + s(phyto,k=4) + s(topog,k=4), family = poisson, offset = log(effort), data = train)
  lm.CM <- gam(CM ~ s(dist,k=3) + s(phyto,k=3) + s(temp_sd,k=4) + s(temp,k=4) + s(salt,k=4), family = poisson, offset = log(effort), data = train)
  lm.GL <- gam(GL ~ s(bathy, k=3) + s(dist,k=3) + s(temp,k=4) + s(phyto,k=4), family = poisson, offset = log(effort), data = train)
  ## apply these models to test data
  test$HP.p <- predict(lm.HP, newdata = test, type = "response", se = FALSE)
  test$HS.p <- predict(lm.HS, newdata = test, type = "response", se = FALSE)
  test$CM.p <- predict(lm.CM, newdata = test, type = "response", se = FALSE)
  test$GL.p <- predict(lm.GL, newdata = test, type = "response", se = FALSE)
  ## compute and store validation metrics
  # prediction error and adj. r-squared
  test$HS[test$HS == 0] <- 0.001 # prevent infinite prediction error
  test$HP[test$HP == 0] <- 0.001 # prevent infinite prediction error
  test$CM[test$CM == 0] <- 0.001 # prevent infinite prediction error
  test$GL[test$GL == 0] <- 0.001 # prevent infinite prediction error
  test <- test %>% mutate(
    d.HP = (((HP - HP.p) / HP) * 100),
    d.HS = (((HS - HS.p) / HS) * 100),
    d.CM = (((CM - CM.p) / CM) * 100),
    d.GL = (((GL - GL.p) / GL) * 100)
  )
  data <- data.frame(
    Zone = z,
    Year = rep(c(2017,2018,2019,2020,2021), times = 4),
    species = rep(c("HP", "HS", "CM", "GL"), each = 5),
    R.squared = rep(c(summary(lm.HP)$r.sq, summary(lm.HS)$r.sq, summary(lm.CM)$r.sq, summary(lm.GL)$r.sq), each = 5),
    Dev.Expl = rep(c(summary(lm.HP)$dev.ex, summary(lm.HS)$dev.ex, summary(lm.CM)$dev.ex, summary(lm.GL)$dev.ex), each = 5)
  )
  LZO_raw <- rbind(LZO_raw, test[,c(1,2,13:20)])
  LZO_results <- rbind(LZO_results, data)
  print(z)
  rm(data, train, test, lm.HP, lm.HS, lm.CM, lm.GL)
}


