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
  filter(Species_code == 'CoMu') %>% 
  mgcv::gam(Count ~  s(bathy, bs = 'bs', m=c(3,1), k = 5) + s(salt, bs = 'bs', m=c(3,1)) + s(dth, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')

## HSeal:
HSeal_daily_beta <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(dist, bs = 'bs', m=c(3,1), k=6) + s(year, bs="re") + s(cruise.gen, bs = "re"),
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
    cruise.gen = rep(1:7, each = 435),
    Effort_Sqkm = rep(1, times = 3045)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

GLFine2_tcur <- 
  cbind(GLFine_tcur,
        predict(GL_daily_beta, newdata = GLFine_tcur, type = "link", se = TRUE, exclude = 's(year)'))

GLFine2_tcur <- 
  GLFine2_tcur %>% 
  dplyr::select(-c(year)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    mod.fit.offset = mod.fit + log(Effort_Sqkm),
    fit = exp(mod.fit.offset),
    LL = exp(mod.fit.offset - (1.96 * se.fit)),
    UL = exp(mod.fit.offset + (1.96 * se.fit)))

GLfine_tcur_plot <- 
  GLFine2_tcur %>%
  unscale('tcur', ., resolution = 'fine') %>%
  mutate(`Cruise #` = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T)) %>% 
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('tcur', ., resolution = 'fine')), aes(x = tcur, y = Density), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = tcur, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = tcur, y = fit)) +
  # scale_color_viridis(discrete = T) +
  # scale_fill_viridis(discrete = T) +
  xlab("Tidal Current Amplitude (m/s)") +
  ylab(expression("Glaucous-winged Gull Density (indiv./km"^2*")")) +
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
    Effort_Sqkm = rep(1, times = 2765),
    year = rep(2017:2021, each = 79) %>% 
      rep(times = 7),
    cruise.gen = rep(1:7, each = 395)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T),
         cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))

GLFine2_dth <- 
  cbind(GLFine_dth,
        predict(GL_daily_beta, newdata = GLFine_dth, type = "link", se = TRUE, exclude = 's(year)'))

GLFine2_dth <- GLFine2_dth %>% 
  dplyr::select(-c(year)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    mod.fit.offset = mod.fit + log(Effort_Sqkm),
    fit = exp(mod.fit.offset),
    LL = exp(mod.fit.offset - (1.96 * se.fit)),
    UL = exp(mod.fit.offset + (1.96 * se.fit)))

GLfine_dth_plot <-
  GLFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'GL') %>% unscale('dth', ., resolution = 'fine')), aes(x = dth, y = Density), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit)) +
  # scale_y_continuous(limits = c(0,38)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  labs(x = "∆ Tide Height (m)") +
  ylab(expression("Glaucous-winged Gull Density (indiv./km"^2*")")) +
  theme_classic() # gulls are more common at intermediate tides

## Common Murre ------------------------------------------------------------

### fine-scale models -----------------------------------------------------
# dist, salt, dth
CMFine_bathy <- 
  data.frame(
    bathy = seq(-1.7,1.4, 0.05) %>% 
      # number of years
      rep(times = 5),
    salt = mean(daily_mbm_grid$salt) %>%
      rep(times = 315),
    dth = mean(daily_mbm_grid$dth) %>%
      rep(times = 315),
    Effort_sqkm = rep(1, times = 315),
    year = rep(2017:2021, each = 63)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_bathy <- 
  cbind(CMFine_bathy,
        predict(CoMu_daily_beta, newdata = CMFine_bathy, type = "link", se.fit = TRUE, exclude = 's(year)'))

CMFine2_bathy <- CMFine2_bathy %>% 
  dplyr::select(-c(year)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    mod.fit.offset = mod.fit + log(Effort_sqkm), 
    fit = exp(mod.fit.offset),
    LL = exp(mod.fit.offset - (1.96 * se.fit)),
    UL = exp(mod.fit.offset + (1.96 * se.fit)) 
  )


CMfine_bathy_plot <-
CMFine2_bathy %>%
  unscale('bathy', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('bathy', ., resolution = 'fine')), aes(x = bathy*-1, y = Density), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = bathy*-1, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy*-1, y = fit)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  # scale_y_continuous(limits = c(0,200)) +
  labs(x = "Water Depth (m)") +
  ylab(expression("Common Murre Density (indiv./km"^2*")")) +
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
    Effort_Sqkm = rep(1, times = 300),
    year = rep(2017:2021, each = 60)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_salt <- 
  cbind(CMFine_salt,
        predict(CoMu_daily_beta, newdata = CMFine_salt, type = "link", se = TRUE, exclude = 's(year)'))

CMFine2_salt <- CMFine2_salt %>% 
  dplyr::select(-c(year)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    mod.fit.offset = mod.fit + log(Effort_Sqkm),
    fit = exp(mod.fit.offset),
    LL = exp(mod.fit.offset - (1.96 * se.fit)),
    UL = exp(mod.fit.offset + (1.96 * se.fit)))

CMfine_salt_plot <-
  CMFine2_salt %>%
  unscale('salt', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('salt', ., resolution = 'fine')), aes(x = salt, y = Density), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = salt, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = salt, y = fit)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  labs(x = "Sea-Surface Salinity (PSU)") +
  ylab(expression("Density (indiv./km"^2*")")) +
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
    Effort_Sqkm = rep(1, times = 395),
    year = rep(2017:2021, each = 79)) %>% 
  mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021), ordered =T))

CMFine2_dth <- 
  cbind(CMFine_dth,
        predict(CoMu_daily_beta, newdata = CMFine_dth, type = "link", se = TRUE, exclude = 's(year)'))

CMFine2_dth <- CMFine2_dth %>% 
  dplyr::select(-c('year')) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    mod.fit.offset = mod.fit + log(Effort_Sqkm),
    fit = exp(mod.fit.offset),
    LL = exp(mod.fit.offset - (1.96 * se.fit)),
    UL = exp(mod.fit.offset + (1.96 * se.fit)))

CMfine_dth_plot <-
  CMFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'CoMu') %>% unscale('dth', ., resolution = 'fine')), aes(x = dth, y = Density), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = dth, y = fit, ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dth, y = fit)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  labs(x = "∆ Tide Height (m)") +
  ylab(expression("Common Murre Density (indiv./km"^2*")")) +
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
        predict(HSeal_daily_beta, newdata = HSFine_dist, type = "link", se = TRUE, exclude = c('s(year)', 's(cruise.gen)')))

HSFine2_dist <- HSFine2_dist %>% 
  dplyr::select(-c(
    # year,
    cruise.gen)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    fit = boot::inv.logit(mod.fit),
    LL = boot::inv.logit(mod.fit - (1.96 * se.fit)),
    UL = boot::inv.logit(mod.fit + (1.96 * se.fit)))

HSfine_dist_plot <-
  HSFine2_dist %>%
  unscale('dist', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
                       mutate(dist = round(dist,0)) %>%
                       group_by(dist, zone, year) %>%
                       summarize(dist = mean(dist),
                                 prob = mean(PresAbs))), aes(x = dist/1000, y = prob), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = dist/1000, y = fit , ymin = LL, ymax = UL), alpha = 0.1) +
  geom_line(aes(x = dist/1000, y = fit)) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_alpha_discrete(range = c(0.5, 1)) +
  # scale_color_manual(values = c("#482878FF", "#5DC863FF")) +
  # scale_fill_manual(values = c("#482878FF", "#5DC863FF")) +
  xlab("Distance from Shore (km)") +
  ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
  guides(color = 'none', fill = 'none') +
  theme_classic() # predicted probability of harbor seal encounter drops to near 0 after 1km from shore

# ggsave(paste0('products/figure3/raw/', Sys.Date(), '_HS_dist.tiff'), device = 'tiff', plot = HSfine_dist_plot, width = 5, height = 4, units = 'in', dpi = 500)

# sst
# HSFine_sst <- 
#   data.frame(
#     sst = seq(-1.9,3.1, 0.1) %>%
#       # number of years
#       rep(times = 5) %>% 
#       # number of cruises
#       rep(times = 7),
#     dist = mean(daily_mbm_grid$dist) %>% rep(times = 1785), 
#     year = rep(2017:2021, each = 51) %>% 
#       rep(times = 7),
#     cruise.gen = rep(1:7, each = 255)) %>% 
#   mutate(year = factor(year, levels = c(2017, 2018, 2019, 2020, 2021)),
#          cruise.gen = factor(cruise.gen, levels = c(1,2,3,4,5,6,7), ordered = T))
# 
# HSFine2_sst <- 
#   cbind(HSFine_sst,
#         predict(HSeal_daily_beta, newdata = HSFine_sst, type = "link", se = TRUE, exclude = c('s(year)', 's(cruise.gen)')))
# 
# HSFine2_sst <- HSFine2_sst %>% 
#   dplyr::select(-c(cruise.gen, year)) %>% 
#   distinct() %>% 
#   mutate(
#     mod.fit = fit,
#     fit = boot::inv.logit(mod.fit),
#     LL = boot::inv.logit(mod.fit - (1.96 * se.fit)),
#     UL = boot::inv.logit(mod.fit + (1.96 * se.fit)))
# 
# HSfine_sst_plot <-
#   HSFine2_sst %>%
#   unscale('sst', ., resolution = 'fine') %>%
#   ggplot() + 
#   # plot effect for a low abundance cruise
#   geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>%
#                        unscale('sst', ., resolution = 'fine') %>%
#                        group_by(zone, cruise.gen) %>%
#                        summarize(sst = mean(sst),
#                                  prob = mean(PresAbs))), aes(x = sst, y = prob), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
#   geom_ribbon(aes(x = sst, y = fit , ymin = LL, ymax = UL), alpha = 0.1) + 
#   geom_line(aes(x = sst, y = fit)) +
#   # scale_x_continuous(limits = c(9.5,10.3)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   # scale_alpha_discrete(range = c(0.5, 1)) +
#   # scale_color_manual(values = c("#482878FF", "#5DC863FF")) +
#   # scale_fill_manual(values = c("#482878FF", "#5DC863FF")) +
#   xlab("Sea Surface Temperature (ºC))") +
#   ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
#   guides(color = 'none', fill = 'none') +
#   theme_classic() # gulls are more common at higher current speeds

# bathy (was highly significant but removed for concurvity)

HSeal_bathy_test <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(bathy, bs = 'bs', m=c(3,1), k=5),
            data = .,
            family = 'binomial')

HSFine_bathy <- 
  data.frame(
    bathy = seq(-1.6,1.4, 0.025))

HSFine2_bathy <- 
  cbind(HSFine_bathy,
        predict(HSeal_bathy_test, newdata = HSFine_bathy, type = "link", se = TRUE))

HSFine2_bathy <- 
  HSFine2_bathy %>% 
  mutate(
    mod.fit = fit,
    fit = boot::inv.logit(mod.fit),
    LL = boot::inv.logit(mod.fit - (1.96 * se.fit)),
    UL = boot::inv.logit(mod.fit + (1.96 * se.fit)))

HSfine_bathy_plot <-
  HSFine2_bathy %>%
  unscale('bathy', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('bathy', ., resolution = 'fine') %>%
                       mutate(dist = round(bathy,0)) %>%
                       group_by(bathy, year) %>%
                       summarize(bahty = mean(bathy),
                                 prob = mean(PresAbs))), aes(x = bathy * -1, y = prob), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = bathy * -1, y = fit , ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = bathy * -1, y = fit)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab("Water Depth (m)") +
  ylab("Predicted Probability of Harbor Seal Encounter (% Chance)") +
  guides(color = 'none', fill = 'none') +
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
        predict(HPorp_daily_beta, newdata = HPFine_dist, type = "link", se = TRUE, exclude = c('s(year)', 's(cruise.gen)')))

HPFine2_dist <- HPFine2_dist %>% 
  dplyr::select(-c(year, cruise.gen)) %>% 
  distinct() %>% 
  mutate(
    mod.fit = fit,
    fit = boot::inv.logit(mod.fit),
    LL = boot::inv.logit(mod.fit - (1.96 * se.fit)),
    UL = boot::inv.logit(mod.fit + (1.96 * se.fit)))

HPfine_dist_plot <-
  HPFine2_dist %>%
  unscale('dist', ., resolution = 'fine') %>%
  ggplot() + 
  # plot effect for a low abundance cruise
  geom_point(data = (daily_mbm_grid %>% filter(Species_code == 'HSeal') %>% unscale('dist', ., resolution = 'fine') %>%
                                            mutate(dist = round(dist,0)) %>%
                                            group_by(dist, year) %>%
                                            summarize(dist = mean(dist),
                                                      prob = mean(PresAbs))), aes(x = dist/1000, y = prob), shape = 21, color = 'black', fill = 'black', alpha = 0.25) +
  geom_ribbon(aes(x = dist/1000, y = fit , ymin = LL, ymax = UL), alpha = 0.1) + 
  geom_line(aes(x = dist/1000, y = fit)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_alpha_discrete(range = c(0.5, 1)) +
  # scale_color_manual(values = c("#482878FF", "#5DC863FF")) +
  # scale_fill_manual(values = c("#482878FF", "#5DC863FF")) +
  # facet_wrap(~Year) +
  xlab("Distance from Shore (km)") +
  ylab("Predicted Probability of Harbor Porpoise Encounter (% Chance)") +
  guides(color = 'none', fill = 'none') +
  theme_classic() # gulls are more common at higher current speeds


# Save Fig.3 --------------------------------------------------------------

zone_comp_plots <-
  list(
    CMfine_bathy_plot,
    CMfine_salt_plot,
    CMfine_dth_plot,
    GLfine_tcur_plot,
    GLfine_dth_plot,
    HSfine_dist_plot,
    # HSfine_sst_plot,
    HSfine_bathy_plot,
    HPfine_dist_plot) %>% 
  set_names(
    c(
      'CM_bathy',
      'CM_salt',
      'CM_dth',
      'GL_tcur',
      'GL_dth',
      'HS_dist',
      # 'HS_sst',
      'HS_bathy',
      'HP_dist'))

for (x in names(zone_comp_plots)) {
  fname <-
    paste0('products/figure4/raw/', Sys.Date(), '_', (x), '.tiff')
  ggsave(fname,
         zone_comp_plots[[x]],
         device = 'tiff',
         width = 6,
         height = 4.5,
         dpi = 500,
         units = 'in')
  print(paste('Done with', x))}

