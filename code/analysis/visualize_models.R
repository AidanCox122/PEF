## THIS CHUNK IS FOR VISUALIZING THE PREDICTOR EFFECTS IN SDMs

# setup -------------------------------------------------------------------

library(tidyverse)
library(PerformanceAnalytics)
library(mgcv)
library(lme4)
library(sf)

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

# create a tibble with scaling information
scaled_env <- 
  env_grid %>% 
  # add delta-tide-height variable
  left_join(dth, by = 'Date') %>% 
  # remove grid cells that do not fall witin zones (cannot be used to train)
  filter(!is.na(zone)) %>% 
  # change zone to a factor
  mutate(zone = factor(zone)) %>% 
  # select variables for model
  dplyr::select(Date, zone, bathy:salt, dth) %>% 
  # find average value in each zone on each cruise date
  group_by(Date, zone) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>%  # n = 1644
  ungroup() %>% 
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
  gam(formula = countInt~ s(bathy,k=3)+s(dist,k=3)+s(phyto,k=3)+s(temp_sd,k=3)+s(sst,k=3),
      family = poisson,
      offset = log(Effort_sqkm),
      data = (interannual_mbm_grid %>% filter(Species_code == 'GL')))
## CM
CoMu_interann_mod <- 
  gam(formula = countInt~s(dist,k=3)+s(bathy,k=3)+s(salt,k=3)+s(phyto,k=3)+s(temp_sd,k=3)+s(sst,k=3),
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
  gam(formula = countInt~ s(sst, k=3)+s(bathy,k=3),
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
  gam(Density ~ s(bathy,k=3)+s(dist, k=3)+s(sst,k=3),
      data = (daily_mbm_grid %>% filter(Species_code == 'GL')),
      family = nb)

## CM:
CoMu_daily_mod <- 
  gam(Density ~ s(dist,k=3)+s(dth, k=3)+s(bathy,k=3)+s(phyto,k=3),
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


# Individual Effect Curves ------------------------------------------------

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

GLcoarse_bathy_plot <- 
  GLCoarse2_bathy %>% 
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('bathy', .)), aes(x = bathy, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Bathymetry (m)") +
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

GLcoarse_dist_plot <- 
  GLCoarse2_dist %>% 
  unscale('dist', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('dist', .)), aes(x = dist/1000, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Distance from Shore (km)") +
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
      sst = seq(-0.75,0.75, 0.025),
      phyto = mean(phyto),
      temp_sd = mean(temp_sd)))

GlCoarse2_sst <- 
  cbind(GlCoarse_sst,
        predict(GL_interann_mod, newdata = GlCoarse_sst, type = "response", se = TRUE))

GlCoarse2_sst <- within(GlCoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GL_coarse_sst_plot <- 
  GlCoarse2_sst %>% 
  unscale('sst', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('sst', .)), aes(x = sst, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea-Surface Temperature (ºF)") +
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

GL_coarse_phyto_plot <- 
  GLCoarse2_phyto %>% 
  unscale('phyto', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('phyto', .)), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Phytoplankton Concentration (µmol/L)") +
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

GL_coarse_tempsd_plot <- 
  GLCoarse2_tempsd %>% 
  unscale('temp_sd', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = temp_sd, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'GL') %>%  unscale('temp_sd', .)), aes(x = temp_sd, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("SST Standard Deviation (ºC)") +
  ylab("Glaucous Gull Count") +
  theme_classic() # gulls are more common in regions high temp_sd


### fine-scale --------------------------------------------------------------
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

GLFine2_bathy %>% 
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('bathy', .)), aes(x = bathy, y = Density), alpha = 0.5) +
  xlab("Bathymetry (m)") +
  ylab("Glaucous Gull Density") +
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

GLFine2_dist %>% 
  unscale('dist', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('dist', .)), aes(x = dist/1000, y = Density), alpha = 0.5) +
  xlab("Distance from Shore (km)") +
  ylab("Glaucous Gull Density") +
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

GLFine2_sst %>% 
  unscale('sst', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('sst', .)), aes(x = sst, y = Density, color = factor(zone)), alpha = 0.5) +
  xlab("Sea-Surface Temperature (ºC)") +
  ylab("Glaucous Gull Density") +
  theme_classic() # gulls are more common at lower sea-surface temperature values


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
  unscale('dist', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dist', .)), aes(x = dist/1000, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
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
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('bathy', .)), aes(x = bathy, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Bathymetry (m)") +
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
      salt = seq(-2,1, 0.025),
      sst = mean(sst)))

CMcoarse2_salt <- 
  cbind(CMcoarse_salt,
        predict(CoMu_interann_mod, newdata = CMcoarse_salt, type = "response", se = TRUE))

CMcoarse2_salt <- within(CMcoarse2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_salt_plot <- 
  CMcoarse2_salt %>% 
  unscale('salt', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('salt', .)), aes(x = salt, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Salinity (ppt)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# phytoplankton
CMcoarse_phyto <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = mean(temp_sd),
      phyto = seq(-1,1, 0.025),
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
  unscale('phyto', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('phyto',.)), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Phytoplankton (µmol/L)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values

# sea-surface temperature standard deviaiton
CMcoarse_tempsd <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      dist = mean(dist),
      temp_sd = seq(-1,1, 0.025),
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
  unscale('temp_sd', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = temp_sd, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = temp_sd, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('temp_sd', .)), aes(x = temp_sd, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
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
      sst = seq(-0.75,0.75, 0.025)))

CMcoarse2_sst <- 
  cbind(CMcoarse_sst,
        predict(CoMu_interann_mod, newdata = CMcoarse_sst, type = "response", se = TRUE))

CMcoarse2_sst <- within(CMcoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CM_coarse_sst_plot <- 
  CMcoarse2_sst %>% 
  unscale('sst', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'CoMu') %>%  unscale('sst', .)), aes(x = sst, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea Surface Temperature (ºC)") +
  ylab("Common Murre Count") +
  theme_classic() # gulls are more common at lower sea-surface temperature values


### fine-scale models -----------------------------------------------------
# dist, dth, bathy, phyto

# distance from shore
CMFine_dist <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = seq(-1,2.5, 0.025),
      dth = mean(dth),
      bathy = mean(bathy),
      phyto = mean(phyto)))

CMFine2_dist <- 
  cbind(CMFine_dist,
        predict(CoMu_daily_mod, newdata = CMFine_dist, type = "response", se = TRUE))

CMFine2_dist <- within(CMFine2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMFine2_dist %>% 
  unscale('dist', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = dist/1000, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dist', .)), aes(x = dist/1000, y = Density), alpha = 0.5) +
  xlab("Distance from Shore (km)") +
  ylab("Common Murre Density") +
  theme_classic() # murres are more common farther from shore

# delta-tide height
CMFine_dth <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = mean(dist),
      dth = seq(-3,2, 0.025),
      bathy = mean(bathy),
      phyto = mean(phyto)))

CMFine2_dth <- 
  cbind(CMFine_dth,
        predict(CoMu_daily_mod, newdata = CMFine_dth, type = "response", se = TRUE))

CMFine2_dth <- within(CMFine2_dth, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMFine2_dth %>% 
  unscale('dth', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dth', .)), aes(x = dth, y = Density, color = zone), alpha = 0.5) +
  xlab("Delta Tide Height (m)") +
  ylab("Common Murre Density") +
  theme_classic() # murres are more common at intermediate tides

# bathymetry
CMFine_bathy <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = mean(dist),
      dth = mean(dth),
      bathy = seq(-2,2, 0.025),
      phyto = mean(phyto)))

CMFine2_bathy <- 
  cbind(CMFine_bathy,
        predict(CoMu_daily_mod, newdata = CMFine_bathy, type = "response", se = TRUE))

CMFine2_bathy <- within(CMFine2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMFine2_bathy %>% 
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('bathy', .)), aes(x = bathy, y = Density), alpha = 0.5) +
  xlab("Bathymetry (m)") +
  ylab("Common Murre Density") +
  theme_classic() # murres are more common in deeper water

# phytoplankton
CMFine_phyto <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = mean(dist),
      dth = mean(dth),
      bathy = mean(bathy),
      phyto = seq(-2,4, 0.025)))

CMFine2_phyto <- 
  cbind(CMFine_phyto,
        predict(CoMu_daily_mod, newdata = CMFine_phyto, type = "response", se = TRUE))

CMFine2_phyto <- within(CMFine2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMFine2_phyto %>% 
  unscale('phyto', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('phyto', .)), aes(x = phyto, y = Density), alpha = 0.5) +
  xlab("Phytoplankton Concentration (µmol/L)") +
  ylab("Common Murre Density") +
  theme_classic() # murres are more common in more productive waters


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
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('bathy', .)), aes(x = bathy, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Bathymetry (m)") +
  ylab("Harbor Seal Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

# phytoplankton

HScoarse_phyto <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      phyto = seq(-1,1, 0.025)))

HScoarse2_phyto <- 
  cbind(HScoarse_phyto,
        predict(HSeal_interann_mod, newdata = HScoarse_phyto, type = "response", se = TRUE))

HScoarse2_phyto <- within(HScoarse2_phyto, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HS_coarse_phyto_plot <- 
HScoarse2_phyto %>% 
  unscale('phyto', .) %>%
  ggplot() + 
  geom_ribbon(aes(x = phyto, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = phyto, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('phyto', .)), aes(x = phyto, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Phytoplankton Concentration (µmol/L)") +
  ylab("Harbor Seal Count") +
  theme_classic() # seals are more common in intermediate phytoplankton concentrations


### fine-scale models -------------------------------------------------------

# bathy, dist, dth

# bathymetry
HSFine_bathy <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = mean(dist),
      dth = mean(dth),
      bathy = seq(-2,2, 0.025)))

HSFine2_bathy <- 
  cbind(HSFine_bathy,
        predict(HSeal_daily_mod, newdata = HSFine_bathy, type = "response", se = TRUE))

HSFine2_bathy <- within(HSFine2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(HSFine2_bathy) + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% group_by(year, zone, bathy) %>% summarize(prob = mean(PresAbs))), aes(x = bathy, y = prob, color = factor(zone)), alpha = 0.5) +
  xlab("Scaled Bathymetry") +
  ylab("Probability of Harbor Seal Sighting") +
  theme_classic()

# distance from shore
HSFine_dist <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = seq(-1,2, 0.025),
      dth = mean(dth),
      bathy = mean(bathy)))

HSFine2_dist <- 
  cbind(HSFine_dist,
        predict(HSeal_daily_mod, newdata = HSFine_dist, type = "response", se = TRUE))

HSFine2_dist <- within(HSFine2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(HSFine2_dist) + 
  geom_ribbon(aes(x = dist, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit)) +
  geom_point(data = (daily_mbm_grid %>%
                       filter(Species_code == 'HSeal') %>%
                       group_by(year, zone, dist) %>%
                       summarize(prob = mean(PresAbs))),
             aes(x = dist, y = prob, color = factor(zone)), alpha = 0.5) +
  xlab("Scaled Distance from Shore") +
  ylab("Probability of Harbor Seal Sighting") +
  theme_classic()

# delta tide height
HSFine_dth <- 
  with(
    daily_mbm_grid,
    data.frame(
      dist = mean(dist),
      dth = seq(-3,2, 0.025),
      bathy = mean(bathy)))

HSFine2_dth <- 
  cbind(HSFine_dth,
        predict(HSeal_daily_mod, newdata = HSFine_dth, type = "response", se = TRUE))

HSFine2_dth <- within(HSFine2_dth, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

ggplot(HSFine2_dth) + 
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit)) +
  geom_point(data = (daily_mbm_grid %>%
                       filter(Species_code == 'HSeal') %>%
                       group_by(year, zone) %>%
                       summarize(dth = mean(dth),
                                 prob = mean(PresAbs))),
             aes(x = dth, y = prob, color = factor(zone)), alpha = 0.5) +
  xlab("Scaled Delta Tide Height") +
  ylab("Probability of Harbor Seal Sighting") +
  theme_classic()


# Harbor Porpoise ---------------------------------------------------------


## coarse-scale models -----------------------------------------------------

# sst, bathy
# sst
HPcoarse_sst <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = mean(bathy),
      sst = seq(-0.6,0.6, 0.025)))

HPcoarse2_sst <- 
  cbind(HPcoarse_sst,
        predict(HPorp_interann_mod, newdata = HPcoarse_sst, type = "response", se = TRUE))

HPcoarse2_sst <- within(HPcoarse2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HPcoarse_sst_plot <- 
  HPcoarse2_sst %>% 
  unscale('sst', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = sst, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp') %>% unscale('sst', .)), aes(x = sst, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Sea Surface Temperature (ºC)") +
  ylab("Harbor Porpoise Count") +
  theme_classic() # seals are more common in shallow water up to a threshold

# bathy
HPcoarse_bathy <- 
  with(
    interannual_mbm_grid,
    data.frame(
      Effort_sqkm = mean(Effort_sqkm),
      bathy = seq(-2,2, 0.025),
      sst = mean(sst)))

HPcoarse2_bathy <- 
  cbind(HPcoarse_bathy,
        predict(HPorp_interann_mod, newdata = HPcoarse_bathy, type = "response", se = TRUE))

HPcoarse2_bathy <- within(HPcoarse2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HPcoarse_bathy_plot <- 
  HPcoarse2_bathy %>% 
  unscale('bathy', .) %>% 
  ggplot() + 
  geom_ribbon(aes(x = bathy, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy, y = fit)) +
  geom_point(data = (interannual_mbm_grid %>% filter(Species_code == 'HPorp') %>% unscale('bathy', .)), aes(x = bathy, y = countInt), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Bathymetry (m)") +
  ylab("Harbor Porpoise Count") +
  theme_classic() # seals are more common in shallow water up to a threshold



# Save Figure 4 -----------------------------------------------------------
coarse_plots <- 
  list(
    # glaucous-winged gull
    GLcoarse_bathy_plot,
    GLcoarse_dist_plot,
    GL_coarse_sst_plot,
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
      'GL_coarse_sst_plot',
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
  # scale env_grid
  # change zone to a factor
  mutate(zone = factor(zone),
         Year = lubridate::year(Date)) %>% 
  # scale environmental variables
  mutate_at(c('bathy', 'topog', 'dist', 'tcur', 'phyto', 'sst', 'temp_sd', 'salt'), base::scale) %>% 
  group_by(Year, grid_id) %>% 
  summarize_if(is.numeric, mean, na.rm = T) %>% 
  ungroup()

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
            HSeal.sd = sd(HSeal, na.rm = T)) %>% 
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
  predicted_predators$CoMu.m %>% quantile(1, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.9, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.7, na.rm = T), 
  predicted_predators$CoMu.m %>% quantile(0.3, na.rm = T),  
  predicted_predators$CoMu.m %>% quantile(0.1, na.rm = T), 
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
  scale_fill_cmocean(name = "balance", discrete = TRUE) +
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
  scale_fill_cmocean(name = "balance", discrete = TRUE) + 
  ggtitle("Harbor Seal Predicted Habitat (Autumn)") +
  xlab("Longitude (ºE)") +
  ylab("Latitude (ºN)") +
  coord_sf() +
  theme_minimal() +
  theme(legend.position="bottom")
