## THIS CHUNK IS FOR VISUALIZING THE PREDICTOR EFFECTS IN SDMs

# setup -------------------------------------------------------------------

library(tidyverse)
library(PerformanceAnalytics)
library(mgcv)
library(lme4)

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


# course-scale models -----------------------------------------------------

# get coarse-scale training data
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

# train best model for:
## GL
GL_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+s(dist,k=3)+s(sst,k=3)+s(phyto,k=3)+s(temp_sd,k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'GL')))
## CM
CoMu_interann_mod <- 
  gam(formula = countInt~s(dist,k=3)+s(temp_sd,k=3)+s(phyto,k=3)+s(bathy,k=3)+s(salt,k=3)+s(sst,k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu')))

## HSeal
HSeal_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+s(phyto,k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal')))

## HPorp
HPorp_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+s(sst, k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp')))

# fine-scale models -------------------------------------------------------

# get fine-scale training data 
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

# train best fine-scale models for:
## GL:
GL_daily_mod <- 
  gam(Density ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

## CM:
CoMu_daily_mod <- 
  gam(Density ~ s(dist,k=3)+s(sst, k=3)+s(salt,k=3)+s(dth,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'CoMu')),
      family = nb)

## HSeal:
HSeal_daily_mod <- 
  glm(PresAbs ~ bathy + dist + dth,
      data = (daily_mbm_grid %>% filter(Species_code == 'HSeal')),
      family = 'binomial')

# HPorp
HPorp_daily_mod <- 
  glm(PresAbs ~ dist,
      data = (daily_mbm_grid %>% filter(Species_code == 'HPorp')),
      family = 'binomial')


# Glaucous gull -----------------------------------------------------------

## coarse-scale ------------------------------------------------------------
# bathy, dist, sst, phyto, temp_sd

# bathymetry
newdata1 <- 
  with(
  interannual_mbm_grid,
  data.frame(
    Effort_sqkm = mean(Effort_sqkm),
    bathy = seq(-1,1, 0.025),
    dist = mean(dist),
    sst = mean(sst),
    phyto = mean(phyto),
    temp_sd = mean(temp_sd)))

newdata2 <- 
  cbind(newdata1,
        predict(GL_interann_mod, newdata = newdata1, type = "response", se = TRUE))

newdata2 <- within(newdata2, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

ggplot(newdata2) + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = bathy, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Bathymetry") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of higher surface temperatures

# distance from shore
newdata_dist <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = seq(-1,2, 0.025),
      sst = mean(sst),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

newdata2_dist <- 
  cbind(newdata_dist,
        predict(GL_interann_mod, newdata = newdata1, type = "response", se = TRUE))

newdata2_dist <- within(newdata2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(newdata2_dist) + 
  geom_ribbon(aes(x = dist, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = dist, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Distance from Shore") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of higher surface temperatures

# sea-surface temperature
newdata_dist <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = seq(-1,2, 0.025),
      sst = mean(sst),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

newdata2_dist <- 
  cbind(newdata_dist,
        predict(GL_interann_mod, newdata = newdata1, type = "response", se = TRUE))

newdata2_dist <- within(newdata2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(newdata2_dist) + 
  geom_ribbon(aes(x = dist, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = dist, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Distance from Shore") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of higher surface temperatures

