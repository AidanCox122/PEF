## THIS CHUNK IS FOR VISUALIZING THE PREDICTOR EFFECTS IN SDMs

# setup -------------------------------------------------------------------

library(tidyverse)
library(PerformanceAnalytics)
library(mgcv)
library(lme4)
library(sf)
library(viridis)

source('code/functions.R')

# fine-scale models -------------------------------------------------------

# train best fine-scale models for:
## GL:
GL_daily_beta <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

## CM:
CoMu_daily_beta <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

## HSeal:
HSeal_daily_beta <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=5) + s(sst, bs = 'bs', m=c(3,1)) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')

# HPorp
HPorp_daily_beta <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=5) + s(year, bs="re") + s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')


# Model Individual Effect Curves ------------------------------------------------

# this section creates effect curves for each variable in our species distribution models

## Glaucous gull -----------------------------------------------------------

### fine-scale --------------------------------------------------------------
# tcur
GLFine_tcur <- 
  data.frame(
    tcur = seq(-1,1.6, 0.03) %>%
      # number of years
      rep(times = 5) %>% 
      # number of cruises
      rep(times = 7),
    dth = mean(daily_mbm_grid$dth) %>% rep(times = 3045), 
    year = rep(2017:2021, each = 87) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 435)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

GLFine2_tcur <- 
  cbind(GLFine_tcur,
        predict(GL_daily_beta, newdata = GLFine_tcur, type = "response", se = TRUE))

GLFine2_tcur <- within(GLFine2_tcur, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GLfine_tcur_plot <- 
  GLFine2_tcur %>%
  unscale('tcur', ., resolution = 'fine') %>%
  mutate(`Year` = factor(year),
         `Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = tcur, y = fit, ymin = LL, ymax = UL, fill = `Year`), alpha = 0.1) + 
  geom_line(aes(x = tcur, y = fit, color = `Year`)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('tcur', ., resolution = 'fine')), aes(x = tcur, y = Density), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Tidal Current Amplitude (m/s)") +
  ylab("Glaucous-Winged Gull Density (indiv./km2)") +
  theme_classic() # gulls are more common at higher current speeds

# dth
GLFine_dth <- 
  data.frame(
    tcur = mean(daily_mbm_grid$tcur) %>%
      rep(times = 2765),
    dth = seq(-2.7,2.0, 0.06) %>% 
      # number of years
      rep(times = 5) %>% 
      # number of cruises
      rep(times = 7), 
    year = rep(2017:2021, each = 79) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 395)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

GLFine2_dth <- 
  cbind(GLFine_dth,
        predict(GL_daily_beta, newdata = GLFine_dth, type = "response", se = TRUE))

GLFine2_dth <- within(GLFine2_dth, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

GLfine_dth_plot <- 
  GLFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  rename(`Year` = year) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL, fill = Year), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit, color = `Year`)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('dth', ., resolution = 'fine')), aes(x = dth, y = Density), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  labs(x = "∆ Tide Height (m)") +
  ylab("Glaucous-Winged Gull Density (indiv./km2)") +
  theme_classic() # gulls are more common at intermediate tides

## Common Murre ------------------------------------------------------------

### fine-scale models -----------------------------------------------------
# bathy, salt, dth

# bathy
CMFine_bathy <- 
  data.frame(
    bathy = seq(-1.7,1.4, 0.05) %>% 
      # number of years
      rep(times = 5),
    salt = mean(daily_mbm_grid$tcur) %>%
      rep(times = 315),
    dth = mean(daily_mbm_grid$dth) %>%
      rep(times = 315), 
    year = rep(2017:2021, each = 63)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_bathy <- 
  cbind(CMFine_bathy,
        predict(CoMu_daily_beta, newdata = CMFine_bathy, type = "response", se = TRUE))

CMFine2_bathy <- within(CMFine2_bathy, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMfine_bathy_plot <- 
  CMFine2_bathy %>%
  unscale('bathy', ., resolution = 'fine') %>%
  rename(`Year` = year) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = bathy * -1, y = fit, ymin = LL, ymax = UL, fill = Year), alpha = 0.1) + 
  geom_line(aes(x = bathy * -1, y = fit, color = `Year`)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('bathy', ., resolution = 'fine')), aes(x = bathy * -1, y = Density), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  labs(x = " Depth (m)") +
  ylab("Common Murre Density (indiv./km2)") +
  theme_classic() # murres most common at shallower depths

# salt
CMFine_salt <- 
  data.frame(
    bathy = mean(daily_mbm_grid$bathy) %>%
      rep(times = 300),
    salt = seq(-4.6,1.3, 0.1) %>% 
      # number of years
      rep(times = 5),
    dth = mean(daily_mbm_grid$dth) %>%
      rep(times = 300), 
    year = rep(2017:2021, each = 60)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_salt <- 
  cbind(CMFine_salt,
        predict(CoMu_daily_beta, newdata = CMFine_salt, type = "response", se = TRUE))

CMFine2_salt <- within(CMFine2_salt, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMfine_salt_plot <- 
  CMFine2_salt %>%
  unscale('salt', ., resolution = 'fine') %>%
  rename(`Year` = year) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL, fill = Year), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit, color = `Year`)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('salt', ., resolution = 'fine')), aes(x = salt, y = Density), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  labs(x = "Sea-Surface Salinity (PSU)") +
  ylab("Common Murre Density (indiv./km2)") +
  theme_classic()


# dth
CMFine_dth <- 
  data.frame(
    bathy = mean(daily_mbm_grid$bathy) %>%
      rep(times = 395),
    salt = mean(daily_mbm_grid$tcur) %>%
      rep(times = 395),
    dth = seq(-2.7,2.0, 0.06) %>% 
      # number of years
      rep(times = 5), 
    year = rep(2017:2021, each = 79)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_dth <- 
  cbind(CMFine_dth,
        predict(CoMu_daily_beta, newdata = CMFine_dth, type = "response", se = TRUE))

CMFine2_dth <- within(CMFine2_dth, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

CMfine_dth_plot <- 
  CMFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  rename(`Year` = year) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL, fill = Year), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit, color = `Year`)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dth', ., resolution = 'fine')), aes(x = dth, y = Density), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  labs(x = "∆ Tide Height (m)") +
  ylab("Common Murre Density (indiv./km2)") +
  theme_classic() # murres most common at intermediate tides

# check average dth in study area
daily_mbm_grid %>% 
  filter(Species_code == 'CoMu') %>% 
  unscale('dth',., resolution = 'fine') %>% 
  pull(dth) %>% 
  mean()

## Harbor Seals ------------------------------------------------------------

### fine-scale models -------------------------------------------------------

# dist + sst

# distance form shore
HSFine_dist <- 
  data.frame(
    dist = seq(-1,2.2, 0.05) %>%
      # number of years
      rep(times = 5) %>% 
      # number of cruises
      rep(times = 7),
    sst = mean(daily_mbm_grid$sst) %>% rep(times = 2275), 
    year = rep(2017:2021, each = 65) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 325)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

HSFine2_dist <- 
  cbind(HSFine_dist,
        predict(HSeal_daily_beta, newdata = HSFine_dist, type = "response", se = TRUE))

HSFine2_dist <- within(HSFine2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)})

HSfine_dist_plot <- 
  HSFine2_dist %>%
  unscale('dist', ., resolution = 'fine') %>%
  mutate(`Year` = factor(year),
         `Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  filter(Year == 2017 | Year == 2020) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dist, y = fit , ymin = LL, ymax = UL, fill = `Cruise #`), alpha = 0.1) + 
  geom_line(aes(x = dist, y = fit, color = `Cruise #`)) +
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
                       mutate(dist = round(dist,0)) %>%
                       group_by(dist, zone) %>%
                       summarize(dist = mean(dist),
                                 prob = mean(PresAbs))), aes(x = dist, y = prob), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  xlab("Distance from Shore (km)") +
  ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
  theme_classic() # gulls are more common at higher current speeds


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

