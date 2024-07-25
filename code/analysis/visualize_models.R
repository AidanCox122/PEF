## THIS CHUNK IS FOR VISUALIZING THE PREDICTOR EFFECTS IN SDMs

# setup -------------------------------------------------------------------

library(tidyverse)
library(PerformanceAnalytics)
library(mgcv)
library(lme4)
library(sf)
library(cmocean)

source('code/functions.R')

env_grid <- 
  read_csv('data/clean/env_grid.csv')

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(zone = Zone) %>%  
  # change zone to a factor
  mutate(zone = factor(zone),
         year = lubridate::year(Date)) %>% 
  # select only data from 2017 onwards
  filter(year >= 2017)

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

# create a tibble with coarse scaling information
coarse_scaled_env <- 
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
  # calculate average conditions in each zone in each year
  group_by(year, zone, Species_code) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>% 
  ungroup() %>% 
  mutate(countInt = round(Count,0)) %>% 
  # scale predictors
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'), base::scale)

vars <- c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt')

coarse_scale_factors <- tibble()
for(x in vars) {
  scale <- 
    coarse_scaled_env %>% 
    pull(x) %>% 
    attr(., 'scaled:scale')
  center <-
    coarse_scaled_env %>% 
    pull(x) %>% 
    attr(., 'scaled:center')
  # store the scaling data in a tibble
  scale_info <- 
    tibble(
      variable = c(x),
      scale = scale,
      center = center)
  # join to repository
  coarse_scale_factors <-
    rbind(coarse_scale_factors, scale_info)
  #cleanup
  rm(scale_info)
  print(x)}


# create a tibble with daily scaling information
scaled_env <- 
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
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth'), base::scale)

vars <- c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt', 'dth')

scale_factors <- tibble()
for(x in vars) {
  scale <- 
    scaled_env %>% 
    pull(x) %>% 
    attr(., 'scaled:scale')
  center <-
    scaled_env %>% 
    pull(x) %>% 
    attr(., 'scaled:center')
  # store the scaling data in a tibble
  scale_info <- 
    tibble(
      variable = c(x),
      scale = scale,
      center = center)
  # join to repository
  scale_factors <-
    rbind(scale_factors, scale_info)
  #cleanup
  rm(scale_info)
  print(x)}

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
  # calculate average conditions in each zone in each year
  group_by(year, zone, Species_code) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>% 
  ungroup() %>% 
  mutate(countInt = round(Count,0)) %>% 
  # scale predictors
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'), base::scale)

# train best model for:
## GL
GL_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+s(dist,k=3)+s(salt,k=5)+s(phyto,k=4)+s(temp_sd,k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'GL')))
## CM
CoMu_interann_mod <- 
  gam(formula = countInt~s(dist,k=3)+s(phyto,k=4)+s(temp_sd,k=4)+s(salt,k=4)+s(sst,k=4)+s(bathy,k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu')))

## HSeal
HSeal_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+s(phyto,k=4),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal')))

## HPorp
HPorp_interann_mod <- 
  gam(formula = countInt~ s(bathy,k=3)+ s(sst,k=4) + s(salt, k=4),
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

# remove attributes attached to each variable by scale function
daily_mbm_grid <-
  unclass(lapply(daily_mbm_grid, function(x) { attributes(x) <- NULL; x })) %>%
  # convert back to a tibble
  data.frame() %>% tibble()


# train best fine-scale models for:
## GL:
GL_daily_mod <- 
  gam(Density ~ s(salt,k=3)+s(zone, bs='re')+s(cruise.gen, bs='re'),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

## CM:
CoMu_daily_mod <- 
  gam(Density ~ s(dth, k = 4) + s(phyto, k = 4) + s(zone, bs = "re") + 
        s(cruise.gen, bs = "re"),
      data = (daily_mbm_grid %>% filter(Species_code == 'CoMu')),
      family = nb)

## HSeal:
HSeal_daily_mod <- 
  glmer(PresAbs ~ dth + (1 | zone),
      data = (daily_mbm_grid %>% filter(Species_code == 'HSeal')),
      family = 'binomial')

# HPorp
HPorp_daily_mod <- 
  glm(PresAbs ~ temp_sd,
      data = (daily_mbm_grid %>% filter(Species_code == 'HPorp')),
      family = 'binomial')


# Coarse-Model Individual Effect Curves ------------------------------------------------

# this section creates effect curves for each variable in our species distribution models

## Glaucous gull -----------------------------------------------------------

### coarse-scale ------------------------------------------------------------
# bathy, dist, phyto, temp_sd, sst

# bathymetry
GLCoarse_bathy <- 
  with(
  interannual_mbm_grid,
  data.frame(
    Effort_sqkm = mean(Effort_sqkm),
    bathy = seq(-2,2, 0.025),
    dist = mean(dist),
    salt = mean(salt),
    phyto = mean(phyto),
    temp_sd = mean(temp_sd)))

GLCoarse2_bathy <- 
  cbind(GLCoarse_bathy,
        predict(GL_interann_mod, newdata = GLCoarse_bathy, type = "response", se = TRUE))

GLCoarse2_bathy <- within(GLCoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

GLcoarse_bathy_plot <- 
  GLCoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy * -1, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy * -1, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('bathy', ., resolution = 'coarse')), aes(x = bathy * -1, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Depth (m)") +
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
      salt = mean(salt),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

GLCoarse2_dist <- 
  cbind(GLCoarse_dist,
        predict(GL_interann_mod, newdata = GLCoarse_dist, type = "response", se = TRUE))

GLCoarse2_dist <- within(GLCoarse2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GLcoarse_dist_plot <- 
  GLCoarse2_dist %>% 
  unscale('dist', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('dist', ., resolution = 'coarse')), aes(x = dist/1000, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Distance from Shore (km)") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common farther from shore

# sea-surface temperature
GlCoarse_salt <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      dist = mean(dist),
      salt = seq(-2.5,1.5, 0.025),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

GlCoarse2_salt <- 
  cbind(GlCoarse_salt,
        predict(GL_interann_mod, newdata = GlCoarse_salt, type = "response", se = TRUE))

GlCoarse2_salt <- within(GlCoarse2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GL_coarse_salt_plot <- 
  GlCoarse2_salt %>% 
  unscale('salt', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('salt', ., resolution = 'coarse')), aes(x = salt, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea-Surface Salinity (PPT)") +
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
      salt = mean(salt),
      phyto = seq(-2,1.8, 0.025),
      temp_sd = mean(temp_sd)))

GLCoarse2_phyto <- 
  cbind(GLCoarse_phyto,
        predict(GL_interann_mod, newdata = GLCoarse_phyto, type = "response", se = TRUE))

GLCoarse2_phyto <- within(GLCoarse2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GL_coarse_phyto_plot <- 
  GLCoarse2_phyto %>% 
  unscale('phyto', ., resolution='coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('phyto', ., resolution='coarse')), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Chlorophyll Concentration (µmol/L)") +
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
      salt = mean(salt),
      phyto = mean(phyto),
      temp_sd = seq(-1.7,2.5, 0.025)))

GLCoarse2_tempsd <- 
  cbind(GLCoarse_tempsd,
        predict(GL_interann_mod, newdata = GLCoarse_tempsd, type = "response", se = TRUE))

GLCoarse2_tempsd <- within(GLCoarse2_tempsd, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GL_coarse_tempsd_plot <- 
  GLCoarse2_tempsd %>% 
  unscale('temp_sd', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = temp_sd, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>%  unscale('temp_sd', ., resolution = 'coarse')), aes(x = temp_sd, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("SST Standard Deviation (ºC)") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions high temp_sd


### fine-scale --------------------------------------------------------------
# salt

# salt
GLFine_salt <- 
  data.frame(
    salt = seq(-4.6,1.3, 0.075) %>%
      # number of zones
      rep(times = 6) %>% 
      # number of cruises
      rep(times = 7),
    zone = rep(1:6, each = 79) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 474))

GLFine2_salt <- 
  cbind(GLFine_salt,
        predict(GL_daily_mod, newdata = GLFine_salt, type = "response", se = TRUE))

GLFine2_salt <- within(GLFine2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

GLfine_salt_plot <- 
  GLFine2_salt %>%
  unscale('salt', ., resolution = 'fine') %>%
  mutate(`Zone` = factor(zone)) %>% 
  filter(cruise.gen == 1 | cruise.gen == 6) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = salt, y = countInt, ymin = LL, ymax = UL, fill = `Zone`), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit, color = `Zone`)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('salt', ., resolution = 'fine') %>% filter(cruise.gen == 1 | cruise.gen == 6)), aes(x = salt, y = Count), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  facet_wrap(~cruise.gen) +
  xlab("Salinity (PPT)") +
  ylab("Glaucous Gull Density (indiv./km2)") +
  theme_classic() # gulls are more common in regions of shallow water

## Common Murre ------------------------------------------------------------


### coarse-scale models -----------------------------------------------------
# dist, bathy, salt, phyto, temp_sd, sst

# distance from shore
CMcoarse_dist <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = seq(-1,2, 0.025),
      temp_sd = mean(temp_sd),
      phyto = mean(phyto),
      bathy = mean(bathy),
      salt = mean(salt),
      sst = mean(sst)))

CMcoarse2_dist <- 
  cbind(CMcoarse_dist,
        predict(CoMu_interann_mod, newdata = CMcoarse_dist, type = "response", se = TRUE))

CMcoarse2_dist <- within(CMcoarse2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_dist_plot <- 
  CMcoarse2_dist %>% 
  unscale('dist', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dist', ., resolution = 'coarse')), aes(x = dist/1000, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Distance from Shore (km)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# bathymetry
CMcoarse_bathy <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = mean(temp_sd),
      phyto = mean(phyto),
      bathy = seq(-2,2, 0.025),
      salt = mean(salt),
      sst = mean(sst)))

CMcoarse2_bathy <- 
  cbind(CMcoarse_bathy,
        predict(CoMu_interann_mod, newdata = CMcoarse_bathy, type = "response", se = TRUE))

CMcoarse2_bathy <- within(CMcoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_bathy_plot <- 
  CMcoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy*-1, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy*-1, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('bathy', ., resolution = 'coarse')), aes(x = bathy*-1, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Depth (m)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# salinity
CMcoarse_salt <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = mean(temp_sd),
      phyto = mean(phyto),
      bathy = mean(bathy),
      salt = seq(-2.5,1.31, 0.025),
      sst = mean(sst)))

CMcoarse2_salt <- 
  cbind(CMcoarse_salt,
        predict(CoMu_interann_mod, newdata = CMcoarse_salt, type = "response", se = TRUE))

CMcoarse2_salt <- within(CMcoarse2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_salt_plot <- 
  CMcoarse2_salt %>% 
  unscale('salt', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('salt', ., resolution = 'coarse')), aes(x = salt, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Salinity (ppt)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# find the salinity which maximizes prediction
CMcoarse2_salt %>%
  unscale('salt', ., resolution = 'coarse') %>%
  mutate(salt = round(salt, 1)) %>% group_by(salt) %>%
  summarize(meanFit = mean(fit)) %>%
  View()
# find average salinity across transect
interannual_mbm_grid %>%
  unscale('salt', ., resolution = 'coarse') %>%
  pull(salt) %>%
  mean()

# phytoplankton
CMcoarse_phyto <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = mean(temp_sd),
      phyto = seq(-2,1.76, 0.025),
      bathy = mean(bathy),
      salt = mean(salt),
      sst = mean(sst)))

CMcoarse2_phyto <- 
  cbind(CMcoarse_phyto,
        predict(CoMu_interann_mod, newdata = CMcoarse_phyto, type = "response", se = TRUE))

CMcoarse2_phyto <- within(CMcoarse2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_phyto_plot <- 
  CMcoarse2_phyto %>% 
  unscale('phyto', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('phyto',., resolution = 'coarse')), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Chlorophyll Concentration (µmol/L)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# find the phytoplankton which maximizes prediction
CMcoarse2_phyto %>%
  unscale('phyto', ., resolution = 'coarse') %>%
  mutate(phyto = round(phyto, 2)) %>% group_by(phyto) %>%
  summarize(meanFit = mean(fit)) %>%
  View()
# find average phytoplankton across transect
interannual_mbm_grid %>%
  unscale('phyto', ., resolution = 'coarse') %>%
  pull(phyto) %>%
  mean()

# sea-surface temperature standard deviaiton
CMcoarse_tempsd <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = seq(-1.7,2.5, 0.025),
      phyto = mean(phyto),
      bathy = mean(bathy),
      salt = mean(salt),
      sst = mean(sst)))

CMcoarse2_tempsd <- 
  cbind(CMcoarse_tempsd,
        predict(CoMu_interann_mod, newdata = CMcoarse_tempsd, type = "response", se = TRUE))

CMcoarse2_tempsd <- within(CMcoarse2_tempsd, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_tempsd_plot <- 
  CMcoarse2_tempsd %>% 
  unscale('temp_sd', ., resolution='coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = temp_sd, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('temp_sd', ., resolution = 'coarse')), aes(x = temp_sd, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("SST Standard Deviation (ºC)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# sea-surface temperature 
CMcoarse_sst <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = mean(temp_sd),
      phyto = mean(phyto),
      bathy = mean(bathy),
      salt = mean(salt),
      sst = seq(-1.65,2.35, 0.025)))

CMcoarse2_sst <- 
  cbind(CMcoarse_sst,
        predict(CoMu_interann_mod, newdata = CMcoarse_sst, type = "response", se = TRUE))

CMcoarse2_sst <- within(CMcoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_sst_plot <- 
  CMcoarse2_sst %>% 
  unscale('sst', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>%  unscale('sst', ., resolution = 'coarse')), aes(x = sst, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea Surface Temperature (ºC)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# find the temperature which maximizes prediction
CMcoarse2_sst %>%
  unscale('sst', ., resolution = 'coarse') %>%
  mutate(sst = round(sst, 1)) %>% group_by(sst) %>%
  summarize(meanFit = mean(fit)) %>%
  View()
# find average salinity across transect
interannual_mbm_grid %>%
  unscale('sst', ., resolution = 'coarse') %>%
  pull(sst) %>%
  mean()

# quick summary of average dth in study area
# daily_mbm_grid %>% unscale('dth', .) %>% filter(Species_code == 'CoMu') %>% pull(dth) %>% mean()

### fine-scale models -----------------------------------------------------
# dth, phyto

# delta-tide height
CMFine_dth <- 
  with(
    daily_mbm_grid,
    data.frame(
      dth = seq(-3,2, 0.05) %>% 
        # number of transects
        rep(times = 6) %>% 
        # number of cruises,
        rep(times = 7),
      phyto = mean(phyto),
      zone = rep(1:6, each = 101) %>% 
        # number of cruises
        rep(times = 7),
      cruise.gen = rep(1:7, each = 606)))

CMFine2_dth <- 
  cbind(CMFine_dth,
        predict(CoMu_daily_mod, newdata = CMFine_dth, type = "response", se = TRUE))

CMFine2_dth <- within(CMFine2_dth, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

CMfine_dth_plot <- 
  CMFine2_dth %>% # group_by(cruise.gen) %>% summarize(meanfit = mean(fit)) %>% View()
  unscale('dth', ., resolution = 'fine') %>%
  mutate(`Zone` = factor(zone)) %>% 
  # filter to the lowest abundance and highest abundance cruise 
  filter(cruise.gen == 2 | cruise.gen == 6) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dth, y = countInt, ymin = LL, ymax = UL, fill = `Zone`), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit, color = `Zone`)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dth', ., resolution = 'fine') %>% filter(cruise.gen == 2 | cruise.gen == 6)), aes(x = dth, y = Count), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  facet_wrap(~cruise.gen) +
  xlab("∆ Tide Height (m)") +
  ylab("Common Murre Density (indiv./km2)") +
  theme_classic() # murres are more common at intermediate tides

# check average dth in study area
daily_mbm_grid %>% 
  filter(Species_code == 'CoMu') %>% 
  unscale('dth',., resolution = 'fine') %>% 
  pull(dth) %>% 
  mean()

# phytoplankton
CMFine_phyto <- 
  with(
    daily_mbm_grid,
    data.frame(
      dth = mean(dth),
      phyto = seq(-2,4, 0.05) %>%
        # number of zones
        rep(times = 6) %>% 
        # number of cruises
        rep(times = 7),
      zone = rep(1:6, each = 121) %>% 
        rep(times = 7),
      cruise.gen = rep(1:7, each = 726)))

CMFine2_phyto <- 
  cbind(CMFine_phyto,
        predict(CoMu_daily_mod, newdata = CMFine_phyto, type = "response", se = TRUE))

CMFine2_phyto <- within(CMFine2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(countInt = round(fit, digits = 0))

CMfine_phyto_plot <- 
  CMFine2_phyto %>%  #group_by(cruise.gen) %>% summarize(meanfit = mean(fit)) %>% View()
  unscale('phyto', ., resolution = 'fine') %>%
  mutate(`Zone` = factor(zone)) %>% 
  # filter to the lowest abundance and highest abundance cruise 
  filter(cruise.gen == 1 | cruise.gen == 6) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = phyto, y = countInt, ymin = LL, ymax = UL, fill = `Zone`), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit, color = `Zone`)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('phyto', ., resolution = 'fine') %>% filter(cruise.gen == 1 | cruise.gen == 6)), aes(x = phyto, y = Count), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  facet_wrap(~cruise.gen) +
  xlab("Chlorophyll Concentration (µmol/L)") +
  ylab("Common Murre Density (indiv./km2)") +
  theme_classic()

# quick summary of average dth in study area
daily_mbm_grid %>% unscale('phyto', .) %>% filter(Species_code == 'CoMu') %>% pull(phyto) %>% mean()
daily_mbm_grid %>% unscale('phyto', .) %>% filter(Species_code == 'CoMu') %>% pull(phyto) %>% sd()

## Harbor Seals ------------------------------------------------------------

### coarse-scale models -----------------------------------------------------

# bathy, phyto
# bathymetry

HScoarse_bathy <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = seq(-2,2, 0.025),
      phyto = mean(phyto)))

HScoarse2_bathy <- 
  cbind(HScoarse_bathy,
        predict(HSeal_interann_mod, newdata = HScoarse_bathy, type = "response", se = TRUE))

HScoarse2_bathy <- within(HScoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HS_coarse_bathy_plot <- 
  HScoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy*-1, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy*-1, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('bathy', ., resolution = 'coarse')), aes(x = bathy*-1, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Depth (m)") +
  ylab("Harbor Seal Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

# find bathymetry that maximizes prediction
HScoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  mutate(bathy = round(bathy, 0)) %>% 
  group_by(bathy) %>% 
  summarize(meanFit = mean(fit)) %>% View()

# find average depth of each zone
interannual_mbm_grid %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  group_by(zone) %>% 
  summarise(meanDepth = mean(bathy))

# phytoplankton
HScoarse_phyto <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      phyto = seq(-2,1.76, 0.025)))

HScoarse2_phyto <- 
  cbind(HScoarse_phyto,
        predict(HSeal_interann_mod, newdata = HScoarse_phyto, type = "response", se = TRUE))

HScoarse2_phyto <- within(HScoarse2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HS_coarse_phyto_plot <- 
HScoarse2_phyto %>% 
  unscale('phyto', ., resolution = 'coarse') %>%
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('phyto', ., resolution = 'coarse')), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Chlorophyll Concentration (µmol/L)") +
  ylab("Harbor Seal Count") +
  theme_classic() # seals are more common in intermediate phytoplankton concentrations

# find the phytoplankton which maximizes prediction
HScoarse2_phyto %>%
  unscale('phyto', ., resolution = 'coarse') %>%
  mutate(phyto = round(phyto, 2)) %>% group_by(phyto) %>%
  summarize(meanFit = mean(fit)) %>%
  View()
# find average phytoplankton across transect
interannual_mbm_grid %>%
  unscale('phyto', ., resolution = 'coarse') %>%
  pull(phyto) %>%
  mean()

### fine-scale models -------------------------------------------------------

# delta tide height
HSFine_dth <- 
  with(
    daily_mbm_grid,
    data.frame(
      dth = seq(-2.75,2, 0.05) %>% 
        rep(times = 6),
      zone = rep(1:6, each = 96)))

HSFine_dth$fit <- predict(HSeal_daily_mod, newdata = HSFine_dth, type = "response")

# HSFine2_dth <- within(HSFine2_dth, {
#   LL <- fit - (1.96 * se.fit)
#   UL <- fit + (1.96 * se.fit)})

HSfine_dth_plot <- 
  HSFine_dth %>% 
  unscale('dth',., resolution = 'fine') %>% 
  mutate(Zone = factor(zone)) %>% 
  ggplot() + 
  # geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL, fill = factor(zone)), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit, color = Zone)) +
  geom_point(data = (daily_mbm_grid %>%
                       filter(Species_code == 'HSeal') %>%
                       unscale('dth', ., resolution = 'fine') %>%
                       mutate(dth = round(dth,0)) %>% 
                       group_by(dth, zone) %>%
                       summarize(dth = mean(dth),
                                 prob = mean(PresAbs))),
             aes(x = dth, y = prob), alpha = 0.5) +
  xlab("Delta Tide Height (m)") +
  ylab("Probability of Harbor Seal Sighting") +
  theme_classic()


# Harbor Porpoise ---------------------------------------------------------


## coarse-scale models -----------------------------------------------------

# sst, bathy, salt
# sst
HPcoarse_sst <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      sst = seq(-1.6,2.32, 0.025),
      salt = mean(salt)))

HPcoarse2_sst <- 
  cbind(HPcoarse_sst,
        predict(HPorp_interann_mod, newdata = HPcoarse_sst, type = "response", se = TRUE))

HPcoarse2_sst <- within(HPcoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HPcoarse_sst_plot <- 
  HPcoarse2_sst %>% 
  unscale('sst', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp') %>% unscale('sst', ., resolution = 'coarse')), aes(x = sst, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea Surface Temperature (ºC)") +
  ylab("Harbor Porpoise Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

# SST at highest prediction
HPcoarse2_sst %>% 
  unscale('sst', ., resolution = 'coarse') %>% 
  mutate(sst = round(sst, 2)) %>% 
  group_by(sst) %>% 
  summarize(meanFit = mean(fit)) %>% View()

# mean SST over transect
interannual_mbm_grid %>% unscale('sst', ., resolution = 'coarse') %>%  filter(Species_code == 'HPorp') %>% pull(sst) %>% mean()

# bathy
HPcoarse_bathy <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = seq(-1.7,1.4, 0.025),
      sst = mean(sst),
      salt = mean(salt)))

HPcoarse2_bathy <- 
  cbind(HPcoarse_bathy,
        predict(HPorp_interann_mod, newdata = HPcoarse_bathy, type = "response", se = TRUE))

HPcoarse2_bathy <- within(HPcoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HPcoarse_bathy_plot <- 
  HPcoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy*-1, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy*-1, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp') %>% unscale('bathy', ., resolution = 'coarse')), aes(x = bathy*-1, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Depth (m)") +
  ylab("Harbor Porpoise Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

# Bathymetry at highest prediction
HPcoarse2_bathy %>% 
  unscale('bathy', ., resolution = 'coarse') %>% 
  mutate(bathy = round(bathy, 1)) %>% 
  group_by(bathy) %>% 
  summarize(meanFit = mean(fit)) %>% View()

# salt
HPcoarse_salt <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      sst = mean(sst),
      salt = seq(-2.5,1.3, 0.025)))

HPcoarse2_salt <- 
  cbind(HPcoarse_salt,
        predict(HPorp_interann_mod, newdata = HPcoarse_salt, type = "response", se = TRUE))

HPcoarse2_salt <- within(HPcoarse2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>%
  # transform predictions less than zero to zero
  mutate(LL = if_else(
    LL < 0,
    0,
    LL))

HPcoarse_salt_plot <- 
  HPcoarse2_salt %>% 
  unscale('salt', ., resolution = 'coarse') %>% 
  ggplot() + 
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp') %>% unscale('salt', ., resolution = 'coarse')), aes(x = salt, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea Surface Salinity (PPT)") +
  ylab("Harbor Porpoise Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

### fine-scale models -------------------------------------------------------

# temp_sd
HPFine_tempSD <- 
  data.frame(
      temp_sd = seq(-1.2,4, 0.05))

HPFine_tempSD$fit <- predict(HPorp_daily_mod, newdata = HPFine_tempSD, type = "response")

# HSFine2_dth <- within(HSFine2_dth, {
#   LL <- fit - (1.96 * se.fit)
#   UL <- fit + (1.96 * se.fit)})

HPfine_tempSD_plot <- 
  HPFine_tempSD %>% 
  unscale('temp_sd',., resolution = 'fine') %>% 
  ggplot() + 
  # geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL, fill = factor(zone)), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (daily_mbm_grid %>%
                       filter(Species_code == 'HPorp') %>%
                       unscale('temp_sd', ., resolution = 'fine') %>%
                       mutate(temp_sd = round(temp_sd,2)) %>% 
                       group_by(temp_sd) %>%
                       summarize(temp_sd = mean(temp_sd),
                                 prob = mean(PresAbs))),
             aes(x = temp_sd, y = prob), alpha = 0.5) +
  xlab("SST Standard Deviation (ºC)") +
  ylab("Probability of Harbor Porpoise Sighting") +
  theme_classic()


# Save Figure 3 -----------------------------------------------------------

fine_plots <- 
  list(
    # glaucous-winged gull
    GLfine_salt_plot,
    # common murre
    CMfine_dth_plot,
    CMfine_phyto_plot,
    # harbor seal
    HSfine_dth_plot,
    # harbor porpoise
    HPfine_tempSD_plot) %>% 
  set_names(
    c(# glaucous-winged gull
      'GLfine_salt_plot',
      # common murre
      'CMfine_dth_plot',
      'CMfine_phyto_plot',
      # harbor seal
      'HSfine_dth_plot',
      # harbor porpoise
      'HPfine_tempSD_plot'))

for (x in names(fine_plots)) {
  fname <-
    paste0('products/figure3/raw/', (x), '.jpg')
  ggsave(fname,
         fine_plots[[x]],
         width = 6,
         height = 4,
         units = 'in')
  print(paste('Done with', x))}

# Save Figure 4 -----------------------------------------------------------
coarse_plots <- 
  list(
    # glaucous-winged gull
    GLcoarse_bathy_plot,
    GLcoarse_dist_plot,
    GL_coarse_salt_plot,
    GL_coarse_phyto_plot,
    GL_coarse_tempsd_plot,
    # common murre
    CM_coarse_dist_plot,
    CM_coarse_bathy_plot,
    CM_coarse_salt_plot,
    CM_coarse_phyto_plot,
    CM_coarse_tempsd_plot,
    CM_coarse_sst_plot,
    # harbor seal
    HS_coarse_bathy_plot,
    HS_coarse_phyto_plot,
    # harbor porpoise
    HPcoarse_sst_plot,
    HPcoarse_bathy_plot) %>% 
  set_names(
    c(# glaucous-winged gull
      'GLcoarse_bathy_plot',
      'GLcoarse_dist_plot',
      'GL_coarse_salt_plot',
      'GL_coarse_phyto_plot',
      'GL_coarse_tempsd_plot',
      # common murre
      'CM_coarse_dist_plot',
      'CM_coarse_bathy_plot',
      'CM_coarse_salt_plot',
      'CM_coarse_phyto_plot',
      'CM_coarse_tempsd_plot',
      'CM_coarse_sst_plot',
      # harbor seal
      'HS_coarse_bathy_plot',
      'HS_coarse_phyto_plot',
      # harbor porpoise
      'HPcoarse_sst_plot',
      'HPcoarse_bathy_plot'))

for (x in names(coarse_plots)) {
  fname <-
    paste0('products/figure4/raw/', (x), '.jpg')
  ggsave(fname,
         coarse_plots[[x]],
         width = 3.68,
         height = 2.98,
         units = 'in')
  print(paste('Done with', x))}

# Prediction Spaces -------------------------------------------------------


## scale prediction data ---------------------------------------------------

## make a tibble with daily values of scaled environmental variables across study area
scaled_env_grid <- 
  env_grid %>% 
  # remove outliers in distnace from shore
  # 1.5 * iqr = 3292.59 m
  # upper outlier limit is 6000m from shore
  mutate(dist = if_else(
    dist > 6000,
    6000,
    dist)) %>%
  # remove outliers in temperature standard deviation
  # 1.5 * iqr = 0.04814385 ºC
  # upper outlier limit is 0.4805773ºC 
  mutate(temp_sd = if_else(
    temp_sd > 0.04814385,
    0.04814385,
    dist)) %>%
  # scale env_grid
  # change zone to a factor
  mutate(zone = factor(zone),
         Year = lubridate::year(Date)) %>% 
  group_by(Year, grid_id) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>% 
  ungroup() %>% 
  # scale environmental variables
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'), base::scale)

## make coarse-scale predictions --------------------------------------------------------

# make predictions of species presence abundance across study area

predicted_grid <- 
  scaled_env_grid %>% 
  mutate(
    GL = predict(GL_interann_mod, newdata = scaled_env_grid, type = "response"),
    CoMu = predict(CoMu_interann_mod, newdata = scaled_env_grid, type = "response"),
    HSeal = predict(HSeal_interann_mod, newdata = scaled_env_grid, type = "response"))


# average abundance/presence across all 5-years

predicted_predators <-
  predicted_grid %>% 
  dplyr::select(grid_id, lat_deg, long_deg, GL:HSeal) %>% 
  group_by(grid_id) %>% 
  summarize(lat_deg = mean(lat_deg),
            long_deg = mean(long_deg),
            GL.m = mean(GL, na.rm = T),
            GL.sd = sd(GL, na.rm = T),
            CoMu.m = mean(CoMu, na.rm = T),
            CoMu.sd = sd(CoMu, na.rm = T),
            HSeal.m = mean(HSeal, na.rm = T), 
            HSeal.sd = sd(HSeal, na.rm = T),
            lat = lat_deg,
            lon = long_deg) %>% 
  st_as_sf(coords = c("long_deg", "lat_deg"), crs = 4326)
  
# find the percentile of abundance / presence for each species in each cell
# PERCENTILE #

# GL
GL.per <- c(0,
            predicted_predators$GL.m %>% quantile(0.1, na.rm = T),
            predicted_predators$GL.m %>% quantile(0.3, na.rm = T),
            predicted_predators$GL.m %>% quantile(0.7, na.rm = T),
            predicted_predators$GL.m %>% quantile(0.9, na.rm = T),
            predicted_predators$GL.m %>% quantile(1, na.rm = T))

# CM
CM.per <- c(
  predicted_predators$CoMu.m %>% quantile(0.1, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.3, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.7, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.9, na.rm = T),  
  predicted_predators$CoMu.m %>% quantile(1, na.rm = T), 
  0)

# HS
HS.per <- c(0,
            predicted_predators$HSeal.m %>% quantile(0.1),
            predicted_predators$HSeal.m %>% quantile(0.3),
            predicted_predators$HSeal.m %>% quantile(0.7),
            predicted_predators$HSeal.m %>% quantile(0.9),
            predicted_predators$HSeal.m %>% quantile(1))

predicted_predators <-
  predicted_predators %>%
  mutate(
    GL.p = cut(GL.m, breaks = GL.per, labels = c(-2,-1,0,1,2)), 
    CM.p = cut(CoMu.m, breaks = CM.per, labels = c(-2,-1,0,1,2)),
    HS.p = cut(HSeal.m, breaks = HS.per, labels = c(-2,-1,0,1,2))) %>% 
  # project into NAD83, UTM Zone 10N
  st_transform(predicted_predators, crs = 26910)


## plot results ------------------------------------------------------------

islands <- st_read('data/GIS/land_area.shp')

## GLAUCOUS WINGED GULL
ggplot() +
  geom_sf(data = islands, fill = "grey") +
  geom_sf(data = predicted_predators, aes(fill = GL.p), size = 2.75, shape = 22) +
  scale_fill_cmocean(name = "balance", discrete = TRUE, labels = c("10", "30", "50", "70", "90")) +
  ggtitle("Glaucous W. Gull Predicted Habitat (Autumn)") +
  xlab("Longitude (ºW)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")
## COMMON MURRE ##
ggplot() +
  geom_sf(data = islands, fill = "grey") +
  geom_sf(data = predicted_predators, aes(fill = CM.p), size = 2.75, shape = 22) +
  scale_fill_cmocean(name = "balance", discrete = TRUE, labels = c("10", "30", "50", "70", "90")) +
  ggtitle("Common Murre Predicted Habitat (Autumn)") +
  xlab("Longitude (ºW)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() + 
  theme(legend.position="bottom")
## HARBOR SEAL ##
ggplot() +
  geom_sf(data = islands, fill = "grey") +
  geom_sf(data = predicted_predators, aes(fill = HS.p), size = 2.75, shape = 22) +
  scale_fill_cmocean(name = "balance", discrete = TRUE, labels = c("10", "30", "50", "70", "90")) + 
  ggtitle("Harbor Seal Predicted Habitat (Autumn)") +
  xlab("Longitude (ºE)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")


### all overlap -------------------------------------------------------------

# areas of overlap
predicted_predators %>% 
  mutate(HS.p = as.numeric(HS.p),
         GL.p = as.numeric(GL.p),
         CM.p = as.numeric(CM.p)) %>%  #View()
  # find regions of core habitat for all species
  filter(HS.p >3 &
           CM.p > 3 &
           GL.p > 3) %>% #write_csv('data/clean/all_predator_sparse.csv')
  # create a metric to plot by
  mutate(HS.p = HS.p-3,
         GL.p = GL.p-3,
         CM.p = CM.p-3,
    all.p = mean(c(HS.p, CM.p, GL.p)) %>% round(0)) %>% # View()
  ggplot()+
  geom_sf(data = islands, fill = "grey") +
  geom_sf(aes(fill = factor(all.p)), size = 2.75, shape = 22) +
  scale_fill_cmocean(name = "balance", discrete = TRUE, labels = c("10", "30", "50", "70", "90")) + 
  ggtitle("Predator Predicted Habitat (Autumn)") +
  xlab("Longitude (ºE)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")

### seabird overlap ---------------------------------------------------------
predicted_predators %>% 
  mutate(HS.p = as.numeric(HS.p),
         GL.p = as.numeric(GL.p),
         CM.p = as.numeric(CM.p)) %>%  #View()
  # find regions of core habitat for all species
  filter(CM.p > 3 &
           GL.p > 3) %>% #write_csv('data/clean/all_predator_sparse.csv')
  # create a metric to plot by
  mutate(HS.p = HS.p-3,
         GL.p = GL.p-3,
         CM.p = CM.p-3,
         all.p = mean(c(HS.p, CM.p, GL.p)) %>% round(0)) %>% # View()
  ggplot()+
  geom_sf(data = islands, fill = "grey") +
  geom_sf(aes(fill = factor(all.p)), size = 2.75, shape = 22) +
  scale_fill_cmocean(name = "balance", discrete = TRUE, labels = c("10", "30", "50", "70", "90")) + 
  ggtitle("Predator Predicted Habitat (Autumn)") +
  xlab("Longitude (ºE)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")

### environmental conditions ----

projected_env_grid <- 
  env_grid %>% 
  st_as_sf(coords = c("long_deg", "lat_deg"), crs = 4326) %>% 
  # project into NAD83, UTM Zone 10N
  st_transform(env_grid, crs = 26910)

ggplot(data = projected_env_grid)+
  geom_sf(data = islands, fill = "grey") +
  geom_sf(aes(fill = phyto), size = 2.75, shape = 22) +
  xlab("Longitude (ºE)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")

