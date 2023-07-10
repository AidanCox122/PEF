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
  # find average value in each zone on each cruise date
  group_by(Date, zone) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>%  # n = 1644
  ungroup() %>% 
  right_join(mbm_data, by = c('Date', 'zone')) %>%
  # lots of NAs from years prior to 2017
  filter(!is.na(bathy)) %>% # n = 672
  mutate(year = lubridate::year(Date)) %>% 
  # scale predictors
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'), base::scale) %>%
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
GLCoarse_bathy <- 
  with(
  interannual_mbm_grid,
  data.frame(
    Effort_sqkm = mean(Effort_sqkm),
    bathy = seq(-2,2, 0.025),
    dist = mean(dist),
    sst = mean(sst),
    phyto = mean(phyto),
    temp_sd = mean(temp_sd)))

GLCoarse2_bathy <- 
  cbind(GLCoarse_bathy,
        predict(GL_interann_mod, newdata = GLCoarse_bathy, type = "response", se = TRUE))

GLCoarse2_bathy <- within(GLCoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

ggplot(GLCoarse2_bathy) + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = bathy, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Bathymetry") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of shallow water

# distance from shore
GLCoarse_dist <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = seq(-1,2, 0.025),
      sst = mean(sst),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

GLCoarse2_dist <- 
  cbind(GLCoarse_dist,
        predict(GL_interann_mod, newdata = GLCoarse_dist, type = "response", se = TRUE))

GLCoarse2_dist <- within(GLCoarse2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLCoarse2_dist) + 
  geom_ribbon(aes(x = dist, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = dist, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Distance from Shore") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common farther from shore

# sea-surface temperature
GlCoarse_sst <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = mean(dist),
      sst = seq(-1,1, 0.025),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

GlCoarse2_sst <- 
  cbind(GlCoarse_sst,
        predict(GL_interann_mod, newdata = GlCoarse_sst, type = "response", se = TRUE))

GlCoarse2_sst <- within(GlCoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GlCoarse2_sst) + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = sst, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Sea-Surface Temperature") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of higher surface temperatures

# phytoplankton
GLCoarse_phyto <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = mean(dist),
      sst = mean(sst),
      phyto = seq(-1,1, 0.025),
      temp_sd = mean(temp_sd)))

GLCoarse2_phyto <- 
  cbind(GLCoarse_phyto,
        predict(GL_interann_mod, newdata = GLCoarse_phyto, type = "response", se = TRUE))

GLCoarse2_phyto <- within(GLCoarse2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLCoarse2_phyto) + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = phyto, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled Phytoplankton Concentration") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of low phytoplankton

# temperature standard deviation
GLCoarse_tempsd <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = mean(dist),
      sst = mean(sst),
      phyto = mean(phyto),
      temp_sd = seq(-1,1, 0.025)))

GLCoarse2_tempsd <- 
  cbind(GLCoarse_tempsd,
        predict(GL_interann_mod, newdata = GLCoarse_tempsd, type = "response", se = TRUE))

GLCoarse2_tempsd <- within(GLCoarse2_tempsd, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLCoarse2_tempsd) + 
  geom_ribbon(aes(x = temp_sd, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL')), aes(x = temp_sd, y = countInt, color = zone), alpha = 0.5) +
  xlab("Scaled SST Standard Deviation") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions high temp_sd


## fine-scale --------------------------------------------------------------
# bathy, dist, sst

# bathymetry
GLFine_bathy <- 
  with(
    daily_mbm_grid,
    data.frame(
      bathy = seq(-2,2, 0.025),
      dist = mean(dist),
      sst = mean(sst)))

GLFine2_bathy <- 
  cbind(GLFine_bathy,
        predict(GL_daily_mod, newdata = GLFine_bathy, type = "response", se = TRUE))

GLFine2_bathy <- within(GLFine2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLFine2_bathy) + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL')), aes(x = bathy, y = Density, color = zone), alpha = 0.5) +
  xlab("Scaled Bathymetry") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions of shallower water

# distance from shore
GLFine_dist <- 
  with(
    daily_mbm_grid,
    data.frame(
      bathy = mean(bathy),
      dist = seq(-1,2.5, 0.025),
      sst = mean(sst)))

GLFine2_dist <- 
  cbind(GLFine_dist,
        predict(GL_daily_mod, newdata = GLFine_dist, type = "response", se = TRUE))

GLFine2_dist <- within(GLFine2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLFine2_dist) + 
  geom_ribbon(aes(x = dist, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL')), aes(x = dist, y = Density, color = zone), alpha = 0.5) +
  xlab("Scaled Distance from Shore") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common farther from shore

# sea-surface temperature
GLFine_sst <- 
  with(
    daily_mbm_grid,
    data.frame(
      bathy = mean(bathy),
      dist = mean(dist),
      sst = seq(-2,3, 0.025)))

GLFine2_sst <- 
  cbind(GLFine_sst,
        predict(GL_daily_mod, newdata = GLFine_sst, type = "response", se = TRUE))

GLFine2_sst <- within(GLFine2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(GLFine2_sst) + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL')), aes(x = sst, y = Density, color = zone), alpha = 0.5) +
  xlab("Scaled Sea-Surface Temperature") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values


# Common Murre ------------------------------------------------------------


## coarse-scale models -----------------------------------------------------
# dist, temp_sd, phyto, bathy, salt, sst

# sea-surface temperature
newdata_sst <- 
  with(
    daily_mbm_grid,
    data.frame(
      bathy = mean(bathy),
      dist = mean(dist),
      sst = seq(-2,3, 0.025)))

newcoarse2_sst <- 
  cbind(newdata_sst,
        predict(GL_daily_mod, newdata = newdata_sst, type = "response", se = TRUE))

newcoarse2_sst <- within(newcoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(newcoarse2_sst) + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL')), aes(x = sst, y = Density, color = zone), alpha = 0.5) +
  xlab("Scaled Sea-Surface Temperature") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

