
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
# null is next best

### interspecies effects HSeal ----

get_logit(
  base = c('bathy','dist', 'dth'),
  test = c('CoMu', 'HPorp', 'GL'),
  species = 'All',
  training = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))))
# CoMu improves model

get_logit(
  base = c('bathy','dist', 'dth','CoMu'),
  test = c('HPorp', 'GL'),
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
  glm(PresAbs ~ bathy + dist + dth + CoMu,
      data = (interspeciesComp %>% mutate(PresAbs = if_else(HSeal > 0, 1, 0))),
      family = 'binomial')

summary(HSeal_daily_mod_interspecies)
# 0.238 dev. expl. (232.8 null; 177.38 residual)
# HPorp not a signif. predictor

# steal the best formula
form_HSeal <- as.formula(HSeal_daily_mod$formula)

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
# 0.0319 deviance explained (199.16 null; 192.79 residual)

# w. interspecies
HPorp_daily_mod_interspecies <- 
  glm(PresAbs ~ dist + HSeal,
      data = (interspeciesComp %>% mutate(PresAbs = if_else(HPorp > 0, 1, 0))),
      family = 'binomial')

summary(HPorp_daily_mod_interspecies) # HSeal nearly signif.
# 0.0484 dev. expl. (199.16 null; 189.52 residual)

# steal the best formula
form_HPorp <- as.formula(HPorp_daily_mod$formula)

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
# HPorp is best 

get_gam(
  base = c('bathy', 'dist', 'sst', 'CoMu', 'HPorp'),
  test = c('HSeal'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = GL)))
# null is best 


### best model for GL ----
GL_daily_mod <- 
  gam(Density ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

summary(GL_daily_mod)
#  30% deviance explained

# w. interspecies 
GL_daily_mod_interspecies <- 
  gam(GL ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3)+s(CoMu,k=3),
      data = interspeciesComp,
      family = nb)

summary(GL_daily_mod_interspecies)
#  40.8 deviance explained

# steal the best formula
form_GL <- as.formula(summary(GL_daily_mod)$formula)

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
# dth is the best model

# forward selection 3
get_gam(
  base = c('dist', 'dth'),
  test = c('bathy', 'phyto', 'sst', 'temp_sd', 'salt'),
  species = 'CoMu',
  training = daily_mbm_grid)
# bathymetry is the best model // BEYOND THIS POINT VARIABLES NOT SIGNIF

# forward selection 4
get_gam(
  base = c('dist', 'dth', 'bathy'),
  test = c('phyto', 'sst', 'temp_sd', 'salt'),
  species = 'CoMu',
  training = daily_mbm_grid)
# phytoplankton is the best model

### interspecies effects CoMu ----
get_gam(
  base = c('dist', 'dth', 'bathy', 'phyto'),
  test = c('GL', 'HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
#GL density improves model

get_gam(
  base = c('dist', 'dth', 'bathy', 'phyto', 'GL'),
  test = c('HSeal', 'HPorp'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# HPorp density improves model

get_gam(
  base = c('dist', 'dth', 'bathy', 'phyto', 'GL', 'HPorp'),
  test = c('HSeal'),
  species = 'All',
  training = (interspeciesComp %>% rename('Density' = CoMu)))
# null is best

### best model for CoMu ----
CoMu_daily_mod <- 
  gam(Density ~ s(dist,k=3)+s(dth, k=3)+s(bathy,k=3)+s(phyto,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'CoMu')),
      family = nb)

summary(CoMu_daily_mod)
#  43.3% deviance explained

# steal the best formula
form_CM <- as.formula(summary(CoMu_daily_mod)$formula)

# W. interspecies
CoMu_daily_mod_interspec <- 
  gam(CoMu ~ s(dist,k=3)+s(dth, k=3)+s(bathy,k=3)+s(phyto,k=3)+s(GL,k=3)+s(HPorp,k=3),
      data = interspeciesComp,
      family = nb)

summary(CoMu_daily_mod_interspec)
# 55% deviance explained


# monte-carlo validation --------------------------------------------------

# this section of code performs monte:carlo validation on fine-scale models trained above
# we run 500 replications and resample 75% of data for training and 25% for testing
# resampling is stratefied by year and cruise day

# create an index of cruise dates to survey from
allrows <- 1:nrow(cruises_ocean_time_all)

# list the rows corresponding to each year
rows_2017 <- 1:6
rows_2018 <- 7:12
rows_2019 <- 13:17
rows_2020 <- 18:24
rows_2021 <- 25:28

seabird.eval <- data.frame()
seabird.coef <- data.frame()

marmam.eval <- data.frame()
marmam.coef <- data.frame()

t1 <- Sys.time()

for (i in 1:500) {
  # step 1: compile the training and testing data
  set.seed(i)
  strat.2017 <- sample(rows_2017, replace = F, size = 0.75*length(rows_2017))
  strat.2018 <- sample(rows_2018, replace = F, size = 0.75*length(rows_2018))
  strat.2019 <- sample(rows_2019, replace = F, size = 0.75*length(rows_2019))
  strat.2020 <- sample(rows_2020, replace = F, size = 0.75*length(rows_2020))
  strat.2021 <- sample(rows_2021, replace = F, size = 0.75*length(rows_2021))
  
  train_rows <- c(strat.2017, strat.2018, strat.2019, strat.2020, strat.2021)
  rm(strat.2017, strat.2018, strat.2019, strat.2020, strat.2021)
  test_rows <- allrows[-train_rows]
  
  # partition env. and species data
  train <- 
    cruises_ocean_time_all[train_rows,] %>%
    dplyr::select(-c(ocean_time)) %>% 
    inner_join(daily_mbm_grid)
  test <- 
    cruises_ocean_time_all[test_rows,] %>%
    dplyr::select(-c(ocean_time)) %>% 
    inner_join(daily_mbm_grid)

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
  glm.HS <- glm(formula = form_HSeal, data= HSeal_train, family = "binomial")
  glm0 <- update(glm.HS, . ~ 1)
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
  
  coef <- coef(glm.HS)
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

