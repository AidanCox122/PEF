
# setup -------------------------------------------------------------------

library(tidyverse) 
library(PerformanceAnalytics)
library(mgcv)
library(lme4)

source('code/functions.R')


env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone)

dth <- 
  read_csv('data/clean/dth.csv')

cruises_ocean_time_all <- data.frame(
  Date = c("2017-10-03", "2017-10-10", "2017-10-24", "2017-10-31", "2017-11-07", "2017-11-16", "2018-10-02", "2018-10-11", "2018-10-16", "2018-10-23", "2018-10-30", "2018-11-06", "2019-10-02", "2019-10-19", "2019-10-23", "2019-10-30", "2019-11-05", "2020-10-06", "2020-10-07", "2020-10-15", "2020-10-20", "2020-10-26", "2020-11-02", "2020-11-10", "2021-10-07", "2021-10-12", "2021-10-19", "2021-10-26"),
  ocean_time = c(276, 283, 297, 304, 311, 320, 640, 649, 654, 661, 668, 675, 1005, 1022, 1026, 1033, 1034, 1375, 1376, 1384, 1389, 1395, 1402, 1410, 1741, 1746, 1753, 1760),
  cruise = c(seq(1,28))) %>% 
  mutate(Date = lubridate::date(Date),
         year = lubridate::year(Date),
         cruise.gen = 1) %>% 
  group_by(year) %>% 
  mutate(cruise.gen = cumsum(cruise.gen)) %>% 
  ungroup()

# format training data ----------------------------------------------------

# add environmental data to each mbm obs. 
daily_mbm_grid <- 
  env_grid %>% 
  # add delta-tide-height variable
  left_join(dth, by = 'Date') %>% 
  # remove grid cells that do not fall witin zones (cannot be used to train)
  filter(!is.na(zone)) %>% 
  # change zone to a factor
  mutate(zone = factor(zone)) %>% 
  # select variables for model
  dplyr::select(Date, zone, bathy:salt,dth) %>% 
  # scale predictors
  mutate_if(is.numeric, base::scale) %>% 
  # find average value in each zone on each cruise date
  group_by(Date, zone) %>% 
  # take the mean of environmental conditions in each zone on each cruise
  summarize_if(is.numeric, mean, na.rm = T) %>%  # n = 1644
  ungroup() %>% 
  right_join(
    mbm_data %>%
      # change zone to factor in mbm_data to match env_grid
      mutate(zone = factor(zone)),
    by = c('Date', 'zone')) %>%
  # lots of NAs from years prior to 2017
  filter(!is.na(bathy)) %>% # n = 672
  mutate(year = lubridate::year(Date),
         PresAbs = if_else(Density > 0,
                           1,
                           0)) %>% 
  # add cruise number
  left_join(
    (cruises_ocean_time_all %>% 
      dplyr::select(Date, cruise.gen)),
    by = 'Date')

interspeciesComp <- 
  daily_mbm_grid %>% 
  dplyr::select(Date:dth, Species_code, Density) %>% 
  pivot_wider(names_from = Species_code,
              values_from = Density) %>% 
  mutate(Species_code = 'All')

# model construction ------------------------------------------------------

daily_mbm_grid %>% 
  dplyr::select('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth') %>% 
  chart.Correlation()

# On random effect terms:
# In the context of a linear regression model, a random effect refers to a
# grouping factor or a categorical variable that introduces additional variation
# in the response variable, which is not directly modeled as a fixed effect.
# Random effects account for the variability that arises due to unobserved or
# unmeasured factors associated with the grouping levels.


# Harbor seal (logit) -------------------------------------------------------------

# FS:1
get_logit(
  test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HSeal',
  training = daily_mbm_grid)
# bathy is the best model

# FS:2
get_logit(
  base = c('bathy'),
  test = c('topog', 'dist', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HSeal',
  training = daily_mbm_grid)
# distance from shore is best model

# forward selection 3
get_logit(
  base = c('bathy', 'dist'),
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HSeal',
  training = daily_mbm_grid)
# dth is next best

# forward selection 3
get_logit(
  base = c('bathy', 'dist', 'dth'),
  test = c('phyto', 'sst', 'temp_sd', 'salt'),
  species = 'HSeal',
  training = daily_mbm_grid)
# dth is next best

### interspecies effects HSeal ----

get_logit(
  base = c('bathy','dist', 'dth'),
  test = c('CoMu', 'HPorp', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# HSPorp improves model

get_logit(
  base = c('bathy','dist', 'dth','HPorp'),
  test = c('CoMu', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# null is best

### best model for HSeal ----
HSeal_daily_mod <- 
  glm(PresAbs ~ bathy + dist + dth,
      data = (daily_mbm_grid %>% filter(Species_code == 'HSeal')),
      family = 'binomial')

summary(HSeal_daily_mod)
# 0.2380584 deviance explained

# w. interspecies
HSeal_daily_mod_interspecies <- 
  glm(PresAbs ~ bathy + dist + dth + HPorp,
      data = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))),
      family = 'binomial')

summary(HSeal_daily_mod_interspecies)
# 0.248 dev. expl.
# HPorp not a signif. predictor

## HPorp (logit) -----------------------------------------------------------

# FS:1
get_logit(
  test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HPorp',
  training = daily_mbm_grid)
# distance from shore is the best model

# FS:2
get_logit(
  base = c('dist'),
  test = c('bathy', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HPorp',
  training = daily_mbm_grid)
# temp_sd from shore is the best model // BUT NOT SIGNIF. 

# FS:3
get_logit(
  base = c('dist', 'temp_sd'),
  test = c('bathy', 'salt', 'dth'),
  species = 'HPorp',
  training = daily_mbm_grid)
# null the best model

### interspecies effects HPorp ----

get_logit(
  base = c('dist'),
  test = c('CoMu', 'HSeal', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HPorp > 0, 1, 0))))
# HSeal improves model

get_logit(
  base = c('dist', 'HSeal'),
  test = c('CoMu', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HPorp > 0, 1, 0))))
# null is best

### best model for HPorp ----
HPorp_daily_mod <- 
  glm(PresAbs ~ dist,
      data = (daily_mbm_grid %>% filter(Species_code == 'HPorp')),
      family = 'binomial')

summary(HPorp_daily_mod)
# 0.042 deviance explained

# w. interspecies
HPorp_daily_mod_interspecies <- 
  glm(PresAbs ~ dist + HSeal,
      data = (interspeciesComp %>% mutate(PresAbs = if_else(HPorp > 0, 1, 0))),
      family = 'binomial')

summary(HPorp_daily_mod_interspecies) # HSeal nearly signif.
# 0.047

## glaucous-winged gull (gam) ----------------------------------------------------

# forward selection 1
get_gam(
  test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'GL',
  training = daily_mbm_grid)
# bathymetry is the best model

# forward selection 2
get_gam(
  base = c('bathy'),
  test = c('topog', 'dist', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'GL',
  training = daily_mbm_grid)
# distance from shore is the best model

# forward selection 3
get_gam(
  base = c('bathy', 'dist'),
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'GL',
  training = daily_mbm_grid)
# sst is the best model

# forward selection 4
get_gam(
  base = c('bathy', 'dist', 'sst'),
  test = c('salt', 'dth'),
  species = 'GL',
  training = daily_mbm_grid)
# delta-tide height is the best model // only nearly signif

# forward selection 5
get_gam(
  base = c('bathy', 'dist', 'sst', 'dth'),
  test = c('salt'),
  species = 'GL',
  training = daily_mbm_grid)
# null is the best model

### interspecies effects GL ----

get_gam(
  base = c('bathy', 'dist', 'sst'),
  test = c('CoMu', 'HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))
# CoMu improves model

get_gam(
  base = c('bathy', 'dist', 'sst', 'CoMu'),
  test = c('HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))
# null is best 


### best model for GL ----
GL_daily_mod <- 
  gam(Density ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(GL_daily_mod)
#  31% deviance explained

# w. interspecies 
GL_daily_mod_interspecies <- 
  gam(GL ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3)+s(CoMu,k=3),
      data = interspeciesComp,
      family = nb)

summary(GL_daily_mod_interspecies)
#  37% deviance explained

## common murre (gam) ----------------------------------------------------

# forward selection 1
get_gam(
  test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'CoMu',
  training = daily_mbm_grid)
# distance from shore is the best model

# forward selection 2
get_gam(
  base = c('dist'),
  test = c('bathy', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'CoMu',
  training = daily_mbm_grid)
# sst is the best model

# forward selection 3
get_gam(
  base = c('dist', 'sst'),
  test = c('bathy', 'salt', 'dth'),
  species = 'CoMu',
  training = daily_mbm_grid)
# salinity is the best model // BEYOND THIS POINT VARIABLES NOT SIGNIF

# forward selection 4
get_gam(
  base = c('dist', 'sst', 'salt'),
  test = c('bathy', 'dth'),
  species = 'CoMu',
  training = daily_mbm_grid)
# delta tide height is the best model

# forward selection 5
get_gam(
  base = c('dist', 'sst', 'salt', 'dth'),
  test = c('bathy'),
  species = 'CoMu',
  training = daily_mbm_grid)
# bathymetry is the best model

### interspecies effects CoMu ----
get_gam(
  base = c('dist', 'sst', 'salt', 'dth'),
  test = c('GL', 'HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
#GL density improves model

get_gam(
  base = c('dist', 'sst', 'salt', 'dth', 'GL'),
  test = c('HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# HPorp density improves model

get_gam(
  base = c('dist', 'sst', 'salt', 'dth', 'GL', 'HPorp'),
  test = c('HSeal'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# null is best

### best model for CoMu ----
CoMu_daily_mod <- 
  gam(Density ~ s(dist,k=3)+s(sst, k=3)+s(salt,k=3)+s(dth,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(CoMu_daily_mod)
#  25.7% deviance explained

# W. interspecies
CoMu_daily_mod_interspec <- 
  gam(CoMu ~ s(dist,k=3)+s(sst, k=3)+s(salt,k=3)+s(dth,k=3)+s(GL,k=3)+s(HPorp,k=3),
      data = interspeciesComp,
      family = nb)

summary(CoMu_daily_mod_interspec)
# 57% deviance explained


# monte-carlo validation --------------------------------------------------



