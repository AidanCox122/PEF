
# setup -------------------------------------------------------------------

library(tidyverse) 
library(PerformanceAnalytics)
library(mgcv)
library(lme4)
library(pROC)

source('code/functions.R')


env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone) %>%  
  # change zone to a factor
  mutate(zone = factor(zone))

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
  # scale predictors
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'), base::scale) %>%
  mutate(year = lubridate::year(Date),
         PresAbs = if_else(Density > 0,
                           1,
                           0)) %>% 
  # add cruise number
  left_join(
    (cruises_ocean_time_all %>% 
       dplyr::select(Date, cruise.gen)),
    by = 'Date') # %>% 
# # add a year column
# mutate(year = lubridate::year(Date))

interspeciesComp <- 
  daily_mbm_grid %>% 
  dplyr::select(Date:dth, year, Species_code, Density) %>% 
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
# glaucous-winged gull ----------------------------------------------------

# random variable selection (1)
nb.GL00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL01 <- gam(Density ~ s(zone, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL02 <- gam(Density ~ s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL03 <- gam(Density ~ s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))

# calculate AIC weights
vec_AIC <- AIC(nb.GL00, nb.GL01, nb.GL02, nb.GL03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

## zone is the best 

# random variable selection (2)
nb.GL00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL01 <- gam(Density ~ s(zone, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL02 <- gam(Density ~ s(zone, bs="re") + s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL03 <- gam(Density ~ s(zone, bs="re") + s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))

# calculate AIC weights
vec_AIC <- AIC(nb.GL00, nb.GL01, nb.GL02, nb.GL03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

## cruise number is the best 

# random variable selection (3)
nb.GL00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL01 <- gam(Density ~ s(zone, bs="re") + s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))
nb.GL02 <- gam(Density ~ s(zone, bs="re") + s(cruise.gen, bs="re") + s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'GL'))

# calculate AIC weights
vec_AIC <- AIC(nb.GL00, nb.GL01, nb.GL02)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

## both models perform equally, sticking with previous iteration 


# forward selection 1
get_gamm(
    test = c('tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
    random = c('zone', 'cruise.gen'),
    species = 'GL',
    training = daily_mbm_grid)

# salinity is the best variable

# forward selection 2
get_gamm(
  base = c('salt'),
  test = c('sst', 'temp_sd', 'dth'),
  random = c('zone', 'cruise.gen'),
  species = 'GL',
  training = daily_mbm_grid)

# null is the best model 

---
## from a previous round of tests ##
# forward selection 3
get_gamm(
  base = c('bathy', 'dist'),
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  random = c('year'),
  species = 'GL',
  training = daily_mbm_grid)

# sst is the best model

# forward selection 4
get_gamm(
  base = c('bathy', 'dist', 'sst'),
  test = c('salt', 'dth'),
  random = c('year'),
  species = 'GL',
  training = daily_mbm_grid)

# dth is the best model but not significant

### interspecies effects GL ----

get_gamm(
  base = c('bathy', 'dist', 'sst'),
  test = c('CoMu', 'HSeal', 'HPorp'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))

# Common Murre was highly significant

get_gamm(
  base = c('bathy', 'dist', 'sst', 'CoMu'),
  test = c('HSeal', 'HPorp'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))

# Harbor Porpoise is best

get_gamm(
  base = c('bathy', 'dist', 'sst', 'CoMu', 'HPorp'),
  test = c('HSeal'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))

# null is best

GL_test_mod <- 
  gam(Density ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3)+s(year, bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(GL_test_mod)
# 30% of deviance explained

# w. interspecies 
GL_test_mod_interspecies <- 
  gam(GL ~ s(bathy,k=3)+s(dist,k=3)+s(sst,k=3)+s(CoMu,k=3)+s(year,bs='re'),
      data = interspeciesComp,
      family = nb)

summary(GL_test_mod_interspecies)
#  40.8 deviance explained

# steal the best formula
form_GL <- as.formula(summary(GL_test_mod)$formula)



# common murre ------------------------------------------------------------

# random variable selection (1)
nb.CM00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM01 <- gam(Density ~ s(zone, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM02 <- gam(Density ~ s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM03 <- gam(Density ~ s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))

# calculate AIC weights
vec_AIC <- AIC(nb.CM00, nb.CM01, nb.CM02, nb.CM03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

# random variable selection (2)
nb.CM00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM01 <- gam(Density ~ s(zone, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM02 <- gam(Density ~ s(zone, bs="re")+s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM03 <- gam(Density ~ s(zone, bs="re")+s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))

# calculate AIC weights
vec_AIC <- AIC(nb.CM00, nb.CM01, nb.CM02, nb.CM03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw
## remaining models perform equally, sticking with zone as only random

# forward selection 1
get_gamm(
  test = c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  random = c('zone'),
  species = 'CoMu',
  training = daily_mbm_grid)

# dth from shore is the best model
# year is not a significant random variable

# forward selection 2
get_gamm(
  base = c('dth'),
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'tcur'),
  random = c('zone'),
  species = 'CoMu',
  training = daily_mbm_grid)
# phytoplankton is the best model

# forward selection 3
get_gamm(
  base = c('dth', 'phyto'),
  test = c('tcur'),
  random = c('year'),
  species = 'CoMu',
  training = daily_mbm_grid)

# tidal current produces the best model

---
  
# forward selection 4
get_gamm(
  base = c('dist', 'dth', 'bathy'),
  test = c('phyto', 'sst', 'temp_sd', 'salt'),
  random = c('year'),
  species = 'CoMu',
  training = daily_mbm_grid)

# phytoplankton produces the best model

### interspecies effects CoMu ----

get_gamm(
  base = c('dist', 'dth', 'bathy', 'phyto'),
  test = c('GL', 'HSeal', 'HPorp'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
#GL density improves model

get_gamm(
  base = c('dist', 'dth', 'bathy', 'phyto', 'GL'),
  test = c('HSeal', 'HPorp'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# HPorp density improves model

get_gamm(
  base = c('dist', 'dth', 'bathy', 'phyto', 'GL', 'HPorp'),
  test = c('HSeal'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# Null model is best

### best model for CoMu ----
CoMu_test_mod <- 
  gam(Density ~ s(dist,k=3)+s(dth, k=3)+s(bathy,k=3)+s(year,bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'CoMu')),
      family = nb)

summary(CoMu_test_mod)
#  43.3% deviance explained

# steal the best formula
form_CM <- as.formula(summary(CoMu_test_mod)$formula)

# W. interspecies
CoMu_daily_mod_interspec <- 
  gam(CoMu ~ s(dist,k=3)+s(dth, k=3)+s(bathy,k=3)+s(phyto,k=3)+s(GL,k=3)+s(HPorp,k=3),
      data = interspeciesComp,
      family = nb)

summary(CoMu_daily_mod_interspec)
# 55% deviance explained
