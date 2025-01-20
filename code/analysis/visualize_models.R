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
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(UL = if_else(UL >= 1, 
                      1, 
                      UL), 
         LL = if_else(LL <= 0, 
                      0, 
                      LL))

HSfine_dist_plot <- 
  HSFine2_dist %>%
  unscale('dist', ., resolution = 'fine') %>%
  mutate(`Year` = factor(year),
         `Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  filter(Year == 2017 | Year == 2020) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dist/1000, y = fit , ymin = LL, ymax = UL, fill = `Cruise #`), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit, color = `Cruise #`)) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
  #                      mutate(dist = round(dist,0)) %>%
  #                      group_by(dist, zone) %>%
  #                      summarize(dist = mean(dist),
  #                                prob = mean(PresAbs))), aes(x = dist, y = prob), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~Year) +
  xlab("Distance from Shore (km)") +
  ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
  theme_classic() # predicted probability of harbor seal encounter drops to near 0 after 1km from shore

# sst
# distance form shore
HSFine_sst <- 
  data.frame(
    sst = seq(-1.9,3.1, 0.1) %>%
      # number of years
      rep(times = 5) %>% 
      # number of cruises
      rep(times = 7),
    dist = mean(daily_mbm_grid$dist) %>% rep(times = 1785), 
    year = rep(2017:2021, each = 51) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 255)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

HSFine2_sst <- 
  cbind(HSFine_sst,
        predict(HSeal_daily_beta, newdata = HSFine_sst, type = "response", se = TRUE))

HSFine2_sst <- within(HSFine2_sst, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(UL = if_else(UL >= 1, 
                      1, 
                      UL), 
         LL = if_else(LL <= 0, 
                      0, 
                      LL))

HSfine_sst_plot <- 
  HSFine2_sst %>%
  unscale('sst', ., resolution = 'fine') %>%
  mutate(`Year` = factor(year),
         `Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  filter(Year == 2017 | Year == 2020) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = sst, y = fit , ymin = LL, ymax = UL, fill = `Cruise #`), alpha = 0.1) + 
  geom_line(aes(x = sst, y = fit, color = `Cruise #`)) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
  #                      mutate(dist = round(dist,0)) %>%
  #                      group_by(dist, zone) %>%
  #                      summarize(dist = mean(dist),
  #                                prob = mean(PresAbs))), aes(x = dist, y = prob), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~Year) +
  xlab("Sea Surface Temperature (ºC))") +
  ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
  theme_classic() # gulls are more common at higher current speeds

# Harbor Porpoise ---------------------------------------------------------

# distance form shore
HPFine_dist <- 
  data.frame(
    dist = seq(-1,2.2, 0.05) %>%
      # number of years
      rep(times = 5) %>% 
      # number of cruises
      rep(times = 7),
    year = rep(2017:2021, each = 65) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 325)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

HPFine2_dist <- 
  cbind(HPFine_dist,
        predict(HPorp_daily_beta, newdata = HPFine_dist, type = "response", se = TRUE))

HPFine2_dist <- within(HPFine2_dist, {
  LL <- fit - (1.96 * se.fit)
  UL <- fit + (1.96 * se.fit)}) %>% 
  mutate(UL = if_else(UL >= 1, 
                      1, 
                      UL), 
         LL = if_else(LL <= 0, 
                      0, 
                      LL))

HPfine_dist_plot <- 
  HPFine2_dist %>%
  unscale('dist', ., resolution = 'fine') %>%
  mutate(`Year` = factor(year),
         `Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  filter(Year == 2017 | Year == 2020) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_ribbon(aes(x = dist/1000, y = fit , ymin = LL, ymax = UL, fill = `Cruise #`), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit, color = `Cruise #`)) +
  # geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
  #                      mutate(dist = round(dist,0)) %>%
  #                      group_by(dist, zone) %>%
  #                      summarize(dist = mean(dist),
  #                                prob = mean(PresAbs))), aes(x = dist, y = prob), shape = 21, color = 'black', fill = NA, alpha = 0.5) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~Year) +
  xlab("Distance from Shore (km)") +
  ylab("Predicted Probability of Harbor Porpoise Encounter (% Chance)") +
  theme_classic() # gulls are more common at higher current speeds


