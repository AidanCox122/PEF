
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
  dplyr::select(Date:dth, year, cruise.gen, Species_code, Density) %>% 
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

### interspecies effects GL ----

get_gamm(
  base = c('salt'),
  test = c('CoMu', 'HSeal', 'HPorp'),
  random = c('zone'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))

# Common Murre was highly significant

get_gamm(
  base = c('salt', 'CoMu'),
  test = c('HSeal', 'HPorp'),
  random = c('zone'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))

# null is best

GL_test_mod <- 
  gam(Density ~ s(salt,k=4)+s(zone, bs='re')+s(cruise.gen, bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(GL_test_mod)
# 33.3% of deviance explained

# w. interspecies 
GL_test_mod_interspecies <- 
  gam(GL ~ s(salt,k=4)+s(CoMu,k=3)+s(zone,bs='re')+s(cruise.gen,bs='re'),
      data = interspeciesComp,
      family = nb)

summary(GL_test_mod_interspecies)
#  37.6 deviance explained

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

#cleanup
rm(nb.CM00, nb.CM01, nb.CM02, nb.CM03, vec_AIC, dAIC, AICw)

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
## cruise.gen is the next best

#cleanup
rm(nb.CM00, nb.CM01, nb.CM02, nb.CM03, vec_AIC, dAIC, AICw)

# random variable selection (3)
nb.CM00 <- gam(Density ~ 1, family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM01 <- gam(Density ~ s(zone, bs="re")+s(cruise.gen, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))
nb.CM02 <- gam(Density ~ s(zone, bs="re")+s(cruise.gen, bs="re")+s(year, bs="re"), family = nb, data=daily_mbm_grid %>% filter(Species_code == 'CoMu'))

# calculate AIC weights
vec_AIC <- AIC(nb.CM00, nb.CM01, nb.CM02)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw


# forward selection 1
get_gamm(
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  random = c('zone', 'cruise.gen'),
  species = 'CoMu',
  training = daily_mbm_grid)

# dth is the best model
# year is not a significant random variable

# forward selection 2
get_gamm(
  base = c('dth'),
  test = c('phyto', 'sst', 'temp_sd', 'salt'),
  random = c('zone', 'cruise.gen'),
  species = 'CoMu',
  training = daily_mbm_grid)
# phytoplankton is the best model

### interspecies effects CoMu ----

get_gamm(
  base = c('dth', 'phyto'),
  test = c('GL', 'HSeal', 'HPorp'),
  random = c('zone', 'cruise.gen'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
#GL density improves model

get_gamm(
  base = c('dth', 'phyto', 'GL'),
  test = c('HSeal', 'HPorp'),
  random = c('zone', 'cruise.gen'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# HPorp density improves model

get_gamm(
  base = c('dth', 'phyto', 'GL', 'HPorp'),
  test = c('HSeal'),
  random = c('year'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# Null model is best

### best model for CoMu ----
CoMu_test_mod <- 
  gam(Density ~ s(dth,k=3)+s(phyto,k=3)+s(zone,bs='re')+s(cruise.gen,bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'CoMu')),
      family = nb)

summary(CoMu_test_mod)
#  45.7% deviance explained

# steal the best formula
form_CM <- as.formula(summary(CoMu_test_mod)$formula)

# W. interspecies
CoMu_daily_mod_interspec <- 
  gam(CoMu ~ +s(dth, k=3)+s(phyto,k=3)+s(GL,k=3)+s(HPorp,k=3)+s(zone,bs='re')+s(cruise.gen,bs='re'),
      data = interspeciesComp,
      family = nb)

summary(CoMu_daily_mod_interspec)
# 54.2% deviance explained


# Harbor seal (logit) -------------------------------------------------------------

# random variable selection (1)
nb.HS00 <- glm(PresAbs ~ 1, family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS01 <- glmer(PresAbs ~ (1|zone), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS02 <- glmer(PresAbs ~ (1|year), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS03 <- glmer(PresAbs ~ (1|cruise.gen), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))

# calculate AIC weights
vec_AIC <- AIC(nb.HS00, nb.HS01, nb.HS02, nb.HS03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

#cleanup
rm(nb.HS00, nb.HS01, nb.HS02, nb.HS03, vec_AIC, dAIC, AICw)

# random variable selection (2)
nb.HS00 <- glmer(PresAbs ~ (1|zone), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS01 <- glmer(PresAbs ~ (1|zone) + (1|year), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS02 <- glmer(PresAbs ~ (1|zone) + (1|cruise.gen), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))

# calculate AIC weights
vec_AIC <- AIC(nb.HS00, nb.HS01, nb.HS02)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

#cleanup
rm(nb.HS00, nb.HS01, nb.HS02, vec_AIC, dAIC, AICw)

# random variable selection (2)
nb.HS00 <- glmer(PresAbs ~ (1|zone) + (1|cruise.gen), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))
nb.HS01 <- glmer(PresAbs ~ (1|zone) + (1|cruise.gen) + (1|year), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HSeal'))

# calculate AIC weights
vec_AIC <- AIC(nb.HS00, nb.HS01)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

#cleanup
rm(nb.HS00, nb.HS01, vec_AIC, dAIC, AICw)


# FS:1
get_mixed_logit(
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  random = c('zone', 'cruise.gen'),
  species = 'HSeal',
  training = daily_mbm_grid)
# dth is the best model

# FS:2
get_mixed_logit(
  base = c('dth'),
  test = c('phyto', 'sst', 'temp_sd', 'salt'),
  random = c('zone', 'cruise.gen'),
  species = 'HSeal',
  training = daily_mbm_grid)
# phyto is best but not statistically significant

### interspecies effects HSeal ----

get_mixed_logit(
  base = c('dth'),
  test = c('CoMu', 'HPorp', 'GL'),
  random = c('zone', 'cruise.gen'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# HPorp improves model

get_mixed_logit(
  base = c('dth', 'HPorp'),
  test = c('CoMu', 'GL'),
  random = c('zone', 'cruise.gen'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# CoMu is best

get_mixed_logit(
  base = c('dth', 'HPorp', 'CoMu'),
  test = c('GL'),
  random = c('zone', 'cruise.gen'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# GL is best

### best model for HSeal ----
HSeal_daily_mod <- 
  glmer(PresAbs ~ dth + (1|zone) + (1|cruise.gen),
      data = (daily_mbm_grid %>% filter(Species_code == 'HSeal')),
      family = 'binomial')

summary(HSeal_daily_mod)
# 0.468 deviance explained

# w. interspecies
HSeal_daily_mod_interspecies <- 
  glmer(PresAbs ~ dth + HPorp + CoMu + GL + (1|zone) + (1|cruise.gen),
      data = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))),
      family = 'binomial')

summary(HSeal_daily_mod_interspecies)
# 0.0617 dev. expl. (232.8 null; 177.38 residual)

# steal the best formula
form_HSeal <- formula(HSeal_daily_mod$formula)

# Harbor Porpoise (logit) -------------------------------------------------------------

# random variable selection (1)
nb.HP00 <- glm(PresAbs ~ 1, family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HPorp'))
nb.HP01 <- glmer(PresAbs ~ (1|zone), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HPorp'))
nb.HP02 <- glmer(PresAbs ~ (1|year), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HPorp'))
nb.HP03 <- glmer(PresAbs ~ (1|cruise.gen), family = 'binomial', data=daily_mbm_grid %>% filter(Species_code == 'HPorp'))

# calculate AIC weights
vec_AIC <- AIC(nb.HP00, nb.HP01, nb.HP02, nb.HP03)[,2]
vec_AIC
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

#cleanup
rm(nb.HP00, nb.HP01, nb.HP02, nb.HP03, vec_AIC, dAIC, AICw)

# FS:1
get_logit(
  test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
  species = 'HPorp',
  training = daily_mbm_grid)
# temperature Standard Deviation is the best model

# FS:2
get_logit(
  base = c('temp_sd'),
  test = c('salt'),
  species = 'HPorp',
  training = daily_mbm_grid)
# null is the best

### interspecies effects HSeal ----

get_logit(
  base = c('temp_sd'),
  test = c('CoMu', 'HSeal', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# HSeal improves model but models did not converge

get_logit(
  base = c('temp_sd', 'HSeal'),
  test = c('CoMu', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# null is best

### best model for HSeal ----
HPorp_daily_mod <- 
  glm(PresAbs ~ temp_sd,
        data = (daily_mbm_grid %>% filter(Species_code == 'HPorp')),
        family = 'binomial')

summary(HPorp_daily_mod)
# 0.468 deviance explained

# w. interspecies
HSeal_daily_mod_interspecies <- 
  glmer(PresAbs ~ dth + HPorp + CoMu + GL + (1|zone) + (1|cruise.gen),
        data = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))),
        family = 'binomial')

summary(HSeal_daily_mod_interspecies)
# 0.013 dev. expl. (199.16 null; 196.43 residual)

# steal the best formula
form_HPorp <- formula(HPorp_daily_mod$formula)


