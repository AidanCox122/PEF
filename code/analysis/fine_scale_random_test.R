
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

#cleanup 
rm(nb.GL00, nb.GL01, nb.GL02, nb.GL03, vec_AIC, dAIC, AICw)

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

#cleanup 
rm(nb.GL00, nb.GL01, nb.GL02, nb.GL03, vec_AIC, dAIC, AICw)

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

#cleanup 
rm(nb.GL00, nb.GL01, nb.GL02, vec_AIC, dAIC, AICw)

# forward selection 1
get_gamm(
    test = c('phyto', 'sst', 'temp_sd', 'salt', 'dth'),
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


## best model for GL -------------------------------------------------------

GL_test_mod <- 
  gam(Density ~ s(salt,k=3)+s(zone, bs='re')+s(cruise.gen, bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(GL_test_mod)
# 33.3% of deviance explained

# w. interspecies 
GL_test_mod_interspecies <- 
  gam(GL ~ s(salt,k=3)+s(CoMu,k=3)+s(zone,bs='re')+s(cruise.gen,bs='re'),
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
# null is best

### best model for HSeal ----
HSeal_daily_mod <- 
  glmer(PresAbs ~ dth + (1|zone) + (1|cruise.gen),
      data = (daily_mbm_grid %>% filter(Species_code == 'HSeal')),
      family = 'binomial')

summary(HSeal_daily_mod)
# 0.468 deviance explained

# steal the best formula
form_HSeal <- formula(HSeal_daily_mod)

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
form_HPorp <- formula(HPorp_daily_mod)

# monte-carlo validation --------------------------------------------------

# this section of code performs monte:carlo validation on fine-scale models trained above
# we run 500 replications and resample 75% of data for training and 25% for testing
# resampling is stratefied by year and cruise day

# create an index of cruise dates to survey from
allrows <- 1:nrow(cruises_ocean_time_all)

# list the rows corresponding to each cruise #
rows_cruise1 <- c(1,7,13,18,25)
rows_cruise2 <- c(2,8,14,19,26)
rows_cruise3 <- c(3,9,15,20,27)
rows_cruise4 <- c(4,10,16,21,28)
rows_cruise5 <- c(5,11,17,22)
rows_cruise6 <- c(6,12,18,23)
rows_cruise7 <- c(24)

seabird.eval <- data.frame()
seabird.coef <- data.frame()

marmam.eval <- data.frame()
marmam.coef <- data.frame()

t1 <- Sys.time()

for (i in 1:500) {
  # step 1: compile the training and testing data
  set.seed(i)
  strat.c1 <- sample(rows_cruise1, replace = F, size = 0.75*length(rows_cruise1))
  strat.c2 <- sample(rows_cruise2, replace = F, size = 0.75*length(rows_cruise2))
  strat.c3 <- sample(rows_cruise3, replace = F, size = 0.75*length(rows_cruise3))
  strat.c4 <- sample(rows_cruise4, replace = F, size = 0.75*length(rows_cruise4))
  strat.c5 <- sample(rows_cruise5, replace = F, size = 0.75*length(rows_cruise5))
  strat.c6 <- sample(rows_cruise6, replace = F, size = 0.75*length(rows_cruise6))
  strat.c7 <- rows_cruise7
  
  train_rows <- c(strat.c1, strat.c2, strat.c3, strat.c4, strat.c5, strat.c6, strat.c7)
  rm(strat.c1, strat.c2, strat.c3, strat.c4, strat.c5, strat.c6, strat.c7)
  test_rows <- allrows[-train_rows]
  
  # partition env. and species data
  train <- 
    daily_mbm_grid %>% 
    filter(Date %in% 
             (cruises_ocean_time_all[train_rows,] %>%
                dplyr::select(-c(ocean_time)) %>% 
                pull(Date)))
  test <- 
    daily_mbm_grid %>% 
    filter(Date %in% 
             (cruises_ocean_time_all[test_rows,] %>%
                dplyr::select(-c(ocean_time)) %>% 
                pull(Date)))
  
  print("Finished Step 1")
  # step 2: add predicted values NOTE: SELECT THE DESIRED MODEL ON LINES 3327:3330 AND 3335:3338
  GL_train <- train %>% filter(Species_code == "GL")
  GL_test <- test %>% filter(Species_code == "GL")
  nb.GL <- gam(formula = form_GL,
               data = GL_train,
               family = nb)
  GL_test$predicted <- predict(nb.GL, newdata = GL_test, type = "response")
  
  CM_train <- train %>% filter(Species_code == "CoMu")
  CM_test <- test %>% filter(Species_code == "CoMu")
  nb.CM <- gam(formula = form_CM,
               data = CM_train,
               family = nb)
  CM_test$predicted <- predict(nb.CM, newdata = CM_test, type = "response")
  
  HPorp_train <- train %>% filter(Species_code == "HPorp")
  HPorp_test <- test %>% filter(Species_code == "HPorp")
  glm.HP <- glm(formula = form_HPorp,
                data = HPorp_train,
                family = 'binomial')
  glm0.0 <- update(glm.HP, . ~ 1)
  HPorp_test$predicted <- predict(glm.HP, newdata = HPorp_test, type = "response")
  
  HSeal_train <- train %>% filter(Species_code == "HSeal")
  HSeal_test <- test %>% filter(Species_code == "HSeal")
  glm.HS <- glmer(formula = form_HSeal, data= HSeal_train, family = "binomial")
  glm0 <- update(glm.HS, . ~ 1 + (1 | zone) + (1 | cruise.gen))
  HSeal_test$predicted <- predict(glm.HS, newdata = HSeal_test, type = "response")
  print("Finished Step 2")
  # step 3: assess the accuracy of the models
  # abundance models
  ## calculate summary statistics
  GL.mi <- min(((GL_test$Density - GL_test$predicted)/GL_test$Density) * 100, na.rm = TRUE)
  GL.me <- median(((GL_test$Density - GL_test$predicted)/GL_test$Density) * 100, na.rm = TRUE)
  GL.ma <- max(((GL_test$Density - GL_test$predicted)/GL_test$Density) * 100, na.rm = TRUE)
  GL.rsq <- summary(nb.GL)$r.sq
  GL.dex <- summary(nb.GL)$dev.expl
  
  CM.mi <- min(((CM_test$Density - CM_test$predicted)/CM_test$Density) * 100, na.rm = TRUE)
  CM.me <- median(((CM_test$Density - CM_test$predicted)/CM_test$Density) * 100, na.rm = TRUE)
  CM.ma <- max(((CM_test$Density - CM_test$predicted)/CM_test$Density) * 100, na.rm = TRUE)
  CM.rsq <- summary(nb.CM)$r.sq
  CM.dex <- summary(nb.CM)$dev.expl
  
  data.eval <- tibble(
    Trial = c(i,i), 
    Species_code = c("GL", "CoMu"), 
    min.pe = c(GL.mi, CM.mi), 
    med.pe = c(GL.me, CM.me), 
    max.pe = c(GL.ma, CM.ma), 
    r.squ = c(GL.rsq, CM.rsq), 
    dev.ex = c(GL.dex, CM.dex),
    test = c(nrow(GL_test), nrow(CM_test)))
  
  coef <- nb.GL$coefficients
  data.cov <- stack(coef)
  data.cov$iter <- i
  data.cov$Species_code <- "GL"
  coef <- nb.CM$coefficients
  coef <- stack(coef)
  coef$iter <- i
  coef$Species_code <- "CM"
  data.cov <- rbind(data.cov, coef)
  
  seabird.eval <- rbind(seabird.eval, data.eval)
  seabird.coef <- rbind(seabird.coef, data.cov)
  rm(GL.mi, CM.mi, GL.me, CM.me, GL.ma, CM.ma, GL.rsq, CM.rsq, GL.dex, CM.dex, coef, data.cov)
  print("Finished Step 3a")
  
  # presence absence models
  # calculate area under receiver opperator curve
  HS.t <- ci.auc(HSeal_test$PresAbs, HSeal_test$predicted)
  HP.t <- ci.auc(HPorp_test$PresAbs, HPorp_test$predicted)
  
  # calculate log liklihood
  HS.rsq <- 1-logLik(glm.HS)/logLik(glm0)
  HP.rsq <- 1-logLik(glm.HP)/logLik(glm0.0)
  
  # calculate percentage of deviance explained
  HS.dex <- 1 - (deviance(glm.HS)/deviance(glm0))
  HP.dex <- 1 - (summary(glm.HP)$deviance/ summary(glm.HP)$null.deviance)
  
  data.eval <- data.frame(
    Trial = c(i,i),
    AUC = c(HS.t[2], HP.t[2]),
    L.CI = c(HS.t[1], HP.t[1]),
    H.CI = c(HS.t[3], HP.t[3]),
    r.sq = c(HS.rsq, HP.rsq),
    dev.ex = c(HS.dex, HP.dex),
    Species_code = c("HSeal", "HPorp")
    
  )
  
  coef <- fixef(glm.HS)
  data.cov <- stack(coef)
  data.cov$iter <- i
  data.cov$Species_code <- "HS"
  coef <- glm.HP$coefficients
  coef <- stack(coef)
  coef$iter <- i
  coef$Species_code <- "HP"
  data.cov <- rbind(data.cov, coef)
  
  marmam.eval <- rbind(marmam.eval, data.eval)
  marmam.coef <- rbind(marmam.coef, data.cov)
  rm(HS.t, HP.t, HS.rsq, HP.rsq, HS.dex, HP.dex, data.eval, data.cov, coef)
  print("Finished Step 3b")
  t2 <- Sys.time()
  print(i)
}

t2 <- Sys.time()

t2 - t1

# evaluate CV results for seabirds
seabird.eval %>% 
  group_by(Species_code) %>% 
  summarize(
    Avg.P.Err = mean(med.pe, na.rm = T),
    SD.P.Err = sd(med.pe, na.rm = T),
    Avg.R2 = mean(r.squ, na.rm = T),
    Avg.DvEx = mean(dev.ex, na.rm = T),
    SD.DevEx = sd(dev.ex, na.rm = T)) %>% View()

# evaluate CV results for mammals
marmam.eval %>% 
  group_by(Species_code) %>% 
  summarize(
    Avg.AUC = mean(AUC, na.rm = T),
    SD.AUC = sd(AUC, na.rm = T),
    Avg.R2 = mean(r.sq, na.rm = T),
    Avg.DvEx = mean(dev.ex, na.rm = T),
    SD.DevEx = sd(dev.ex, na.rm = T)) %>% View()



