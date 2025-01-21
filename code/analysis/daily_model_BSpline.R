
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
    by = 'Date') %>% 
  mutate(
    year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021, ordered = T)),
    cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

interspeciesComp <- 
  daily_mbm_grid %>% 
  dplyr::select(Date:dth, year, Species_code, Density) %>% 
  pivot_wider(names_from = Species_code,
              values_from = Density) %>% 
  mutate(Species_code = 'All')

# model construction ------------------------------------------------------

## HSeal -------------------------------------------------------------------
# B-spline
# perform variable selection
daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(topog, bs = 'bs', m=c(3,1), k = 5) + s(dist, bs = 'bs', m=c(3,1), k = 5) + s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(phyto, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(temp_sd, bs = 'bs', m=c(3,1)) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial',
            select = TRUE) %>% 
  summary() # dist, tcur, bathy, topog, sst, and dth

# parameterize the full model
# remove tcur, dth, and topog, not signif
HSeal_daily_beta <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=5) + s(sst, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')


summary(HSeal_daily_beta)

# test the interaction term
HSeal.interact.01 <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=5) + s(sst, bs = 'bs', m=c(3,1)) + s(dist, by = sst, bs = 'bs', m = c(3,1), k=5) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')

summary(HSeal.interact.01) # interaction not significant - ignoring

# test for concurvity
mgcv::concurvity(HSeal_daily_beta, full = T)
# no high concurvity between fixed effects

rm(HSeal_interact_01)

## HPorp -------------------------------------------------------------------
# B- Spline

# perform variable selection
daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(topog, bs = 'bs', m=c(3,1), k = 5) + s(dist, bs = 'bs', m=c(3,1), k = 5) + s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(phyto, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(temp_sd, bs = 'bs', m=c(3,1)) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial',
            select = TRUE) %>% 
  summary() # dist, maybe salt

HPorp_daily_beta <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=5) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')

summary(HPorp_daily_beta)

# assess random effect terms
mgcv::gam.vcomp(HPorp_daily_beta)
mgcv::vcov.gam(HPorp_daily_beta) %>% summary() %>% View()

# test for concurvity
mgcv::concurvity(HPorp_daily_beta, full = T)

## GL ----------------------------------------------------------------------
# B-Spline

# select model terms
daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(topog, bs = 'bs', m=c(3,1), k = 5) + s(dist, bs = 'bs', m=c(3,1), k = 5) + s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(phyto, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(temp_sd, bs = 'bs', m=c(3,1)) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb',
            select = TRUE) %>% 
  summary() # terms are bathy, dist, dth, topog, sst

# define the full model
GL_daily_beta <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(GL_daily_beta)
# 44.8% dev. explained

# test the interaction term
GL.interact.01 <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(dth, bs = 'bs', m=c(3,1)) + s(tcur, by = dth, bs = 'bs', m=c(3,1), k = 5) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(GL.interact.01) # interaction not significant - ignoring

# assess random effect terms
mgcv::gam.vcomp(GL_daily_beta)
mgcv::vcov.gam(GL_daily_beta) %>% summary() %>% View()

# test for concurvity
mgcv::concurvity(GL_daily_beta, full = T)
# very high values all around
mgcv::concurvity(GL_daily_beta, full = F)
# bathy highly concuve with dist, tcur
# dist also highly concurve with tcur

rm(GL_interact_01)

## CoMu --------------------------------------------------------------------
# B-Spline

# select model terms
daily_mbm_grid %>% 
  filter(Species_code == 'CoMu') %>% 
  mgcv::gam(Count ~ s(bathy, bs = 'bs', m=c(3,1), k= 5) + s(topog, bs = 'bs', m=c(3,1), k = 5) + s(dist, bs = 'bs', m=c(3,1), k = 5) + s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(phyto, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(temp_sd, bs = 'bs', m=c(3,1)) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb',
            select = TRUE) %>% 
  summary() # terms are dist, phyto, bathy, salt, dth, topog

# define the full model, pulled out phyto and topog bc. not signif
CoMu_daily_beta <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(CoMu_daily_beta)
# 44.2% dev. explained

# assess random effect terms
mgcv::gam.vcomp(CoMu_daily_beta)
mgcv::vcov.gam(CoMu_daily_beta) %>% View()


# test for concurvity
mgcv::concurvity(CoMu_daily_beta, full = T)
# high concurvity across the board
mgcv::concurvity(CoMu_daily_beta, full = F)
# bathymetry and distance from shore are highly colinear
# removing dist has a smaller negative impact on deviance explained

# test the interaction term
CoMu.interact.01 <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(bathy, by = salt, bs = 'bs', m=c(3,1), k = 5) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(CoMu.interact.01) # interaction significant

CoMu.interact.02 <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(bathy, by = dth, bs = 'bs', m=c(3,1), k = 5) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(CoMu.interact.02) # interaction not significant

CoMu.interact.03 <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(salt, by = dth, bs = 'bs', m=c(3,1), k = 5) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(CoMu.interact.03) # interaction not significant


vec_AIC <- c(AIC(CoMu_daily_beta), AIC(CoMu.interact.01), AIC(CoMu.interact.02), AIC(CoMu.interact.03))
dAIC <- vec_AIC - min(vec_AIC)
AICw <- exp(-dAIC/2) / sum(exp(-dAIC/2))
AICw

rm(CoMu.interact.01, CoMu.interact.02, CoMu.interact.03)

# the results suggest that interaction 3 is best but still not significant

# Tweedie Distribution ----------------------------------------------------

# select model terms
daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(topog, bs = 'bs', m=c(3,1), k = 5) + s(dist, bs = 'bs', m=c(3,1), k = 5) + s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(phyto, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(temp_sd, bs = 'bs', m=c(3,1)) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            offset = log(Effort_sqkm),
            family = Tweedie(p = 1.9),
            select = TRUE) %>% 
  summary() # terms are bathy, dist, dth, topog, sst, tcur

# define the full model
GL_daily_tweedie <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = Tweedie(p = 1.9))

summary(GL_daily_tweedie)
# 44.8% dev. explained

# test for concurvity
mgcv::concurvity(GL_daily_tweedie, full = T)
# high concurvity across the board
mgcv::concurvity(GL_daily_tweedie, full = F)



# Leave-One-Out Validation ------------------------------------------------

## Year --------------------------------------------------------------------

LYO_metrics <- data.frame()
LYO_raw <- data.frame()
for (y in unique(daily_mbm_grid$year)) {
  ## create train and test sets
  counter = 1
  train <- daily_mbm_grid %>% filter(year != y)
  test <- daily_mbm_grid %>% filter(year == y)
  ## train a model for each species
  
  cv_models <- 
    list(
    HPorp = gam(formula = formula.gam(HPorp_daily_beta),
                family = 'binomial',
                offset = log(Effort_sqkm),
                data = train %>% filter(Species_code == 'HPorp')),
    HSeal = gam(formula = formula.gam(HSeal_daily_beta),
                family = 'binomial',
                offset = log(Effort_sqkm),
                data = train %>% filter(Species_code == 'HSeal')),
    CoMu = gam(formula = formula.gam(CoMu_daily_beta),
               family = 'nb',
               offset = log(Effort_sqkm),
               data = train %>% filter(Species_code == 'CoMu')),
   GL = gam(formula = formula.gam(GL_daily_beta),
            family = 'nb',
            offset = log(Effort_sqkm),
            data = train %>% filter(Species_code == 'GL')))
  
  ## apply these models to test data
  
  test.full <- tibble()
  performance.full <- tibble()
  
  for(x in unique(daily_mbm_grid$Species_code)) {
    # filter out data for testing
    test.species <- 
      test %>% filter(Species_code == x) 
    
    # select model of interest
    model <- 
      cv_models[[x]]
    
    # make predictions
    test.predicted <- 
      test.species %>% 
      cbind(
        mgcv::predict.gam(model, newdata = test.species, exclude = "s(year)", type = 'response', se = T)) %>% 
      mutate(
        Pred.Response = fit %>% round(digits = 2),
        # convert 0 predicted response to 0.01 so that Obs:Pred can be calculated
        Pred.Response = if_else(Pred.Response == 0,
                                0.01,
                                Pred.Response),
        Obs.Response = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          PresAbs,
          Density),
        # categorize seabird abundance into high(1) or low(0) based on time-series avg. density 
        Obs.Cat = case_when(
          # assign high low abundance to glaucous gull
          Species_code == 'GL' &
            Obs.Response >= (daily_mbm_grid %>%
                               # this bit pulls out the time-series mean
                               filter(Species_code == 'GL') %>%
                               pull(Density) %>%
                               median()) ~ 1,
          Species_code == 'GL' &
            Obs.Response < (daily_mbm_grid %>%
                              filter(Species_code == 'GL') %>%
                              pull(Density) %>%
                              median()) ~ 0,
          # do the same thing for common murre
          Species_code == 'CoMu' &
            Obs.Response >= (daily_mbm_grid %>%
                               filter(Species_code == 'CoMu') %>%
                               pull(Density) %>%
                               median()) ~ 1,
          Species_code == 'CoMu' &
            Obs.Response < (daily_mbm_grid %>%
                              filter(Species_code == 'CoMu') %>%
                              pull(Density) %>%
                              median()) ~ 0,
          Species_code %in% c("HSeal", "HPorp") ~ NA),
        PresAbs,
        Response.Type = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          'PresAbs',
          'Density')) %>% 
      dplyr::select(-c(Count, Density, PresAbs, fit))
    
    # calculate the proportion of observations marine mammals were present in zone x for each year
    Obs.Prop.Table <- 
      test.predicted %>%
      group_by(Species_code, zone, year) %>%
      summarize(
        Obs.Prop = mean(Obs.Response) %>% round(digits = 2))
    
    test.predicted <- 
      test.predicted %>% 
      left_join(Obs.Prop.Table, by = c('Species_code', 'year', 'zone')) %>% 
      mutate(
        `Obs:Pred` = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          Obs.Prop/Pred.Response,
          Obs.Response/Pred.Response))
    
    
    # gather information for spearman's ranked correlation test
    # rank each zone in order of abundance for spearman's correlation
    
    test.ranks <-
      test.predicted %>% 
      # calculate prop. marine mammal sightings in each zone and average pred. prob. of occurance in each zone
      group_by(Species_code, zone) %>% 
      mutate(
        Obs.Prop = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          mean(Obs.Response) %>% round(digits = 2),
          NA),
        Avg.Pred.Response = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          mean(Pred.Response) %>% round(digits = 2),
          NA)) %>%  #filter(Species_code == 'HSeal') %>% filter(zone == 1 | zone == 6) %>% filter(Obs.Prop == Obs.Prop) %>%  View()
      group_by(Date, Species_code) %>% 
      reframe(
        zone,
        Response.Type,
        Obs.Response,
        Obs.Prop,
        `Obs.Rank` = if_else(
          Species_code %in% c('GL', 'CoMu'),
          rank(Obs.Response, ties.method = 'average'),
          rank(Obs.Prop, ties.method = 'average')),
        Pred.Response,
        Avg.Pred.Response,
        `Pred.Rank` = if_else(
          Species_code %in% c('GL', 'CoMu'),
          rank(Pred.Response, ties.method = 'average'),
          rank(Avg.Pred.Response, ties.method = 'average'))) %>% 
      ungroup()

    
    if(x %in% c('HSeal', 'HPorp')) {
      O <- test.ranks %>%
        filter(Species_code == x) %>%
        pull(Obs.Rank) %>% 
        head()} else{
          O <- test.ranks %>%
            filter(Species_code == x) %>%
            pull(Obs.Rank)}
    
    if(x %in% c('HSeal', 'HPorp')) {
      R <- test.ranks %>% 
        filter(Species_code == x) %>% 
        pull(Pred.Rank) %>% 
        head()} else{
          R <- test.ranks %>% 
            filter(Species_code == x) %>% 
            pull(Pred.Rank)}
    
    test.full <- rbind(test.full, test.predicted)
    
    # calculate the ROC and TSS
    if(x %in% c('GL', 'CoMu')) {
      test.roc <- 
        roc(test.predicted$Obs.Cat, test.predicted$Pred.Response)
      test.auc <- auc(test.predicted$Obs.Cat, test.predicted$Pred.Response)
      } else{
        test.roc <- 
          roc(test.predicted$Obs.Response, test.predicted$Pred.Response)
        test.auc <- auc(test.predicted$Obs.Response, test.predicted$Pred.Response)}
    
    # plot(test.roc)
    
    roc_data <- 
      data.frame(sensitivity = test.roc$sensitivities, specificity = test.roc$specificities, thresholds = test.roc$thresholds) %>% 
      mutate(TSS = (sensitivity + specificity) - 1)
    
    # store the success metrics
    data <- tibble(
      year = y,
      species = x,
      Dev.Expl = summary(model)$dev.expl,
      AUC = test.auc,
      TSS = max(roc_data$TSS),
      `Obs:Pred` = mean(test.predicted$`Obs:Pred`),
      `sd.Obs:Pred` = sd(test.predicted$`Obs:Pred`),
      spearman.rho = cor.test(O, R, method = 'spearman')$estimate,
      spearman.p = cor.test(O, R, method = 'spearman')$p.value)
    
    # save to repo
    performance.full <- 
      rbind(performance.full, data)
  # tell the user you're done with species x
    print(
      paste('Done with', x, sep = ' '))
  }
  
  LYO_raw <- rbind(LYO_raw, test.full[,-c(3:13)])
  LYO_metrics <- rbind(LYO_metrics, performance.full)
  print(y)
  rm(data, train, test, cv_models, performance.full, test.full, test.ranks, test.predicted, y, x, O, R)
}


## Zone --------------------------------------------------------------------

# same process but for zones
LZO_metrics <- data.frame()
LZO_raw <- data.frame()
for (y in unique(daily_mbm_grid$zone)) {
  ## create train and test sets
  counter = 1
  train <- daily_mbm_grid %>% filter(zone != y)
  test <- daily_mbm_grid %>% filter(zone == y)
  ## train a model for each species
  
  cv_models <- 
    list(
      HPorp = gam(formula = formula.gam(HPorp_daily_beta),
                  family = 'binomial',
                  offset = log(Effort_sqkm),
                  data = train %>% filter(Species_code == 'HPorp')),
      HSeal = gam(formula = formula.gam(HSeal_daily_beta),
                  family = 'binomial',
                  offset = log(Effort_sqkm),
                  data = train %>% filter(Species_code == 'HSeal')),
      CoMu = gam(formula = formula.gam(CoMu_daily_beta),
                 family = 'nb',
                 offset = log(Effort_sqkm),
                 data = train %>% filter(Species_code == 'CoMu')),
      GL = gam(formula = formula.gam(GL_daily_beta),
               family = 'nb',
               offset = log(Effort_sqkm),
               data = train %>% filter(Species_code == 'GL')))
  
  ## apply these models to test data
  
  test.full <- tibble()
  performance.full <- tibble()
  
  for(x in unique(daily_mbm_grid$Species_code)) {
    # filter out data for testing
    test.species <- 
      test %>% filter(Species_code == x) 
    
    # select model of interest
    model <- 
      cv_models[[x]]
    
    # make predictions
    test.predicted <- 
      test.species %>% 
      cbind(
        mgcv::predict.gam(model, newdata = test.species, type = 'response', se = T)) %>% 
      mutate(
        Pred.Response = fit %>% round(digits = 2),
        # convert 0 predicted response to 0.01 so that Obs:Pred can be calculated
        Pred.Response = if_else(Pred.Response == 0,
                                0.01,
                                Pred.Response),
        Obs.Response = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          PresAbs,
          Density),
        # categorize seabird abundance into high(1) or low(0) based on time-series avg. density 
        Obs.Cat = case_when(
          # assign high low abundance to glaucous gull
          Species_code == 'GL' &
            Obs.Response >= (daily_mbm_grid %>%
                               # this bit pulls out the time-series mean
                               filter(Species_code == 'GL') %>%
                               pull(Density) %>%
                               median()) ~ 1,
          Species_code == 'GL' &
            Obs.Response < (daily_mbm_grid %>%
                              filter(Species_code == 'GL') %>%
                              pull(Density) %>%
                              median()) ~ 0,
          # do the same thing for common murre
          Species_code == 'CoMu' &
            Obs.Response >= (daily_mbm_grid %>%
                               filter(Species_code == 'CoMu') %>%
                               pull(Density) %>%
                               median()) ~ 1,
          Species_code == 'CoMu' &
            Obs.Response < (daily_mbm_grid %>%
                              filter(Species_code == 'CoMu') %>%
                              pull(Density) %>%
                              median()) ~ 0,
          Species_code %in% c("HSeal", "HPorp") ~ NA),
        PresAbs,
        Response.Type = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          'PresAbs',
          'Density')) %>% 
      dplyr::select(-c(Count, Density, PresAbs, fit))
    
    # calculate the proportion of observations marine mammals were present in zone x for each year
    Obs.Prop.Table <- 
      test.predicted %>%
      group_by(Species_code, zone, year) %>%
      summarize(
        Obs.Prop = mean(Obs.Response) %>% round(digits = 2))
    
    test.predicted <- 
      test.predicted %>% 
      left_join(Obs.Prop.Table, by = c('Species_code', 'year', 'zone')) %>% 
      mutate(
        `Obs:Pred` = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          Obs.Prop/Pred.Response,
          Obs.Response/Pred.Response))
    
    # gather information for spearman's ranked correlation test
    # rank each zone in order of abundance for spearman's correlation
    
    test.ranks <-
      test.predicted %>% 
      # calculate prop. marine mammal sightings in each zone and average pred. prob. of occurance in each zone
      group_by(Species_code, year) %>% 
      mutate(
        Obs.Prop = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          mean(Obs.Response) %>% round(digits = 2),
          NA),
        Avg.Pred.Response = if_else(
          Species_code %in% c('HSeal', 'HPorp'),
          mean(Pred.Response) %>% round(digits = 2),
          NA)) %>%  #View()
      ungroup() %>% 
      transmute(
        Date,
        year,
        zone,
        Species_code,
        Response.Type,
        Obs.Response,
        Obs.Prop,
        `Obs.Rank` = if_else(
          Species_code %in% c('GL', 'CoMu'),
          rank(Obs.Response, ties.method = 'average'),
          rank(Obs.Prop, ties.method = 'average')),
        Pred.Response,
        Avg.Pred.Response,
        `Pred.Rank` = if_else(
          Species_code %in% c('GL', 'CoMu'),
          rank(Pred.Response, ties.method = 'average'),
          rank(Avg.Pred.Response, ties.method = 'average'))) %>% 
      ungroup() #%>% View()
    
    
      O <- test.ranks %>%
      filter(Species_code == x) %>%
      pull(Obs.Rank)
    
    R <- test.ranks %>% 
      filter(Species_code == x) %>% 
      pull(Pred.Rank)
    
    test.full <- rbind(test.full, test.predicted)
    
    # calculate the ROC and TSS
    if(x %in% c('GL', 'CoMu')) {
      test.roc <- 
        roc(test.predicted$Obs.Cat, test.predicted$Pred.Response)
      test.auc <- auc(test.predicted$Obs.Cat, test.predicted$Pred.Response)
    } else{
      test.roc <- 
        roc(test.predicted$Obs.Response, test.predicted$Pred.Response)
      test.auc <- auc(test.predicted$Obs.Response, test.predicted$Pred.Response)}
    
    # plot(test.roc)
    
    roc_data <- 
      data.frame(sensitivity = test.roc$sensitivities, specificity = test.roc$specificities, thresholds = test.roc$thresholds) %>% 
      mutate(TSS = (sensitivity + specificity) - 1)
    
    # store the success metrics
    data <- tibble(
      year = y,
      species = x,
      Dev.Expl = summary(model)$dev.expl,
      AUC = test.auc,
      TSS = max(roc_data$TSS),
      `Obs:Pred` = mean(test.predicted$`Obs:Pred`),
      `sd.Obs:Pred` = sd(test.predicted$`Obs:Pred`),
      spearman.rho = cor.test(O, R, method = 'spearman')$estimate,
      spearman.p = cor.test(O, R, method = 'spearman')$p.value)
    
    # save to repo
    performance.full <- 
      rbind(performance.full, data)
    # tell the user you're done with species x
    print(
      paste('Done with', x, sep = ' '))
  }
  
  LZO_raw <- rbind(LZO_raw, test.full[,-c(3:13)])
  LZO_metrics <- rbind(LZO_metrics, performance.full)
  print(y)
  rm(data, train, test, cv_models, performance.full, test.full, test.ranks, test.predicted, y, x, O, R)
}

# table 5
# LYO
LYO_metrics %>% 
  group_by(species) %>% 
  summarise(
    Avg.Dev.Expl = mean(Dev.Expl),
    SD.Dev.Expl = sd(Dev.Expl),
    # spacer
    Avg.AUC = mean(AUC),
    sd.AUC = sd(AUC),
    # spacer
    Avg.TSS = mean(TSS),
    sd.TSS. = sd(TSS))
# LZO
LZO_metrics %>% 
  group_by(species) %>% 
  summarise(
    Avg.Dev.Expl = mean(Dev.Expl),
    SD.Dev.Expl = sd(Dev.Expl),
    # spacer
    Avg.AUC = mean(AUC),
    sd.AUC = sd(AUC),
    # spacer
    Avg.TSS = mean(TSS),
    sd.TSS. = sd(TSS))
