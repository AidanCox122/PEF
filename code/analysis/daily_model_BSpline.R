
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
HSeal_daily <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(bathy, k=3) + s(dist, k=3) + s(sst,k=3) + s(dth, k=3) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial',
            select = TRUE)

summary(HSeal_daily)

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

# test for concurvity
mgcv::concurvity(HSeal_daily_beta, full = T)
# no high concurvity between fixed effects

## HPorp -------------------------------------------------------------------
# Thin-Plate
HPorp_daily <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(dist, k=3) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial',
            select = TRUE)

summary(HPorp_daily)

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

# define the full model
GL_daily <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(bathy, k=3) + s(dist, k=3) + s(sst,k=3) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(GL_daily)

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

# assess random effect terms
mgcv::gam.vcomp(GL_daily_beta)
mgcv::vcov.gam(GL_daily_beta) %>% summary() %>% View()

# test for concurvity
mgcv::concurvity(GL_daily_beta, full = T)
# very high values all around
mgcv::concurvity(GL_daily_beta, full = F)
# bathy highly concuve with dist, tcur
# dist also highly concurve with tcur

## CoMu --------------------------------------------------------------------

CoMu_daily <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(bathy,k=3) + s(sst,k=3) + s(phyto, k=3) + s(salt, k=3) + s(topog,k=3) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

summary(CoMu_daily)

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

LYO_results <- data.frame()
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
  
  for(y in unique(daily_mbm_grid$Species_code)) {
    # filter out data for testing
    test.species <- 
      test %>% filter(Species_code == y) 
    
    # select model of interest
    model <- 
      cv_models[[y]]
    
    # make predictions
    test.predicted <- 
      test.species %>% 
      cbind(
        mgcv::predict.gam(model, newdata = test.species, exclude = c("s(year)"), newdata.guaranteed = T, type = 'response', se = T)) %>% 
      mutate(
        Dev.Expl = summary(model)$dev.expl,
        Pred.Response = fit,
        Obs.Response = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          PresAbs,
          Density),
        Response.Type = if_else(
          Species_code %in% c("HSeal", "HPorp"),
          'PresAbs',
          'Density'),
        `Obs:Pred` = (Obs.Response/Pred.Response)) %>% 
      dplyr::select(-c(Count, Density, PresAbs, fit))
    
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
          rank(Obs.Response, ties.method = 'average'),
          rank(Avg.Pred.Response, ties.method = 'average'))) %>% 
      ungroup()

    
    if(y %in% c('HSeal', 'HPorp')) {
      O <- test.ranks %>%
        filter(Species_code == y) %>%
        pull(Obs.Rank) %>% 
        head()} else{
          O <- test.ranks %>%
            filter(Species_code == y) %>%
            pull(Obs.Rank)}
    
    if(y %in% c('HSeal', 'HPorp')) {
      R <- test.ranks %>% 
        filter(Species_code == y) %>% 
        pull(Pred.Rank) %>% 
        head()} else{
          R <- test.ranks %>% 
            filter(Species_code == y) %>% 
            pull(Pred.Rank)}
    
    # extract rho coefficient
    test.predicted$spearman.rho <- cor.test(O, R, method = 'spearman')$estimate
    # extract p-value
    test.predicted$spearman.p <- cor.test(O, R, method = 'spearman')$p.value
    
    
    test.full <- rbind(test.full, test.predicted)
    
    print(
      paste('Done with', y, sep = ' '))
  }
  
  data <- data.frame(
    year = y,
    species = c("HP", "HS", "CM", "GL"),
    R.squared = c(summary(lm.HP)$r.sq, summary(lm.HS)$r.sq, summary(lm.CM)$r.sq, summary(lm.GL)$r.sq),
    Dev.Expl = c(summary(lm.HP)$dev.ex, summary(lm.HS)$dev.ex, summary(lm.CM)$dev.ex, summary(lm.GL)$dev.ex)
  )
  LYO_raw <- rbind(LYO_raw, test.full[,c(1,2,13:24)])
  LYO_results <- rbind(LYO_results, data)
  print(y)
  rm(data, train, test, lm.HP, lm.HS, lm.CM, lm.GL)
}

# same process but for zones
LZO_results <- data.frame()
LZO_raw <- data.frame()
for (z in unique(cv_base$zone)) {
  ## create train and test sets
  train <- cv_base %>% filter(zone != z)
  test <- cv_base %>% filter(zone == z)
  ## train a model for each species
  lm.HP <- gam(formula = HPorp~s(bathy,k=3)+ s(sst,k=4) + s(salt, k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.HS <- gam(formula = HSeal~ s(bathy,k=3)+s(phyto,k=4),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.CM <- gam(formula = CoMu~s(dist,k=3)+s(phyto,k=4)+s(temp_sd,k=4)+s(salt,k=4)+s(sst,k=4)+s(bathy,k=3),
               family = poisson,
               offset = log(Effort_sqkm),
               data = train)
  lm.GL <- gam(formula = GL~ s(bathy,k=3)+s(dist,k=3)+s(salt,k=5)+s(phyto,k=4)+s(temp_sd,k=4),
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
      # cap prediction at highest observed value
      HP.p = if_else(HP.p > max(train$HPorp),
                     max(train$HPorp),
                     HP.p),
      # repeat for Harbor seal
      HS.p = predict(lm.HS, newdata = test, type = "response", se = FALSE),
      HS.p = round(HS.p, digits = 0),
      HS.p = if_else(HS.p > max(train$HSeal),
                     max(train$HSeal),
                     HS.p),
      # Common Murre
      CM.p = predict(lm.CM, newdata = test, type = "response", se = FALSE),
      CM.p = round(CM.p, digits = 0),
      CM.p = if_else(CM.p > max(train$CoMu),
                     max(train$CoMu),
                     CoMu),
      # Glaucous-winged gull
      GL.p = predict(lm.GL, newdata = test, type = "response", se = FALSE),
      GL.p = round(GL.p, digits = 0),
      GL.p = if_else(GL.p > max(train$GL),
                     max(train$GL),
                     GL.p)) %>% 
    # # convert zero values to very small numbers to avoid infinite predition error
    # mutate_all( ~ ifelse( . == 0, 0.000001, .)) %>%
    # calculate prediction error for each species in each zone
    mutate(
      d.HP = (((HPorp - HP.p) / HPorp) * 100),
      d.HS = (((HSeal - HS.p) / HSeal) * 100),
      d.CM = (((CoMu - CM.p) / CoMu) * 100),
      d.GL = (((GL - GL.p) / GL) * 100))
  
  data <- data.frame(
    zone = z,
    species = c("HP", "HS", "CM", "GL"),
    R.squared = c(summary(lm.HP)$r.sq, summary(lm.HS)$r.sq, summary(lm.CM)$r.sq, summary(lm.GL)$r.sq),
    Dev.Expl = c(summary(lm.HP)$dev.ex, summary(lm.HS)$dev.ex, summary(lm.CM)$dev.ex, summary(lm.GL)$dev.ex)
  )
  
  LZO_raw <- rbind(LZO_raw, test.full[,c(1,2,13:24)])
  LZO_results <- rbind(LZO_results, data)
  print(z)
  rm(data, train, test, lm.HP, lm.HS, lm.CM, lm.GL)
}