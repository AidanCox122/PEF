# exploring covariate effects

# GL ----------------------------------------------------------------------

# tcur
GLFine2_tcur %>%
  unscale('tcur', ., resolution = 'fine') %>%
  filter(tcur <= 4) %>% 
  group_by(year) %>% 
  filter(fit == max(fit)) %>% View()


daily_mbm_grid %>%
  unscale('tcur', ., resolution = 'fine') %>%
  filter(Species_code == 'GL') %>% 
  dplyr::select(zone, tcur) %>% 
  distinct()

# dth
GLFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  filter(dth >= 1.8 & dth <= 2.25) %>%
  group_by(year) %>% 
  filter(fit == max(fit)) %>% View()

daily_mbm_grid %>% 
  unscale('dth', ., resolution = 'fine') %>%
  filter(Species_code == 'GL') %>% 
  group_by(year) %>% 
  summarize(
    mean.dth = mean(dth))
  
# variability explianed by null model for gulls
GL_year <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
  mgcv::gam(Count ~ s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')
summary(GL_year)


# CoMu --------------------------------------------------------------------

# bathy
CMFine2_bathy %>%
  unscale('bathy', ., resolution = 'fine') %>%
  # filter(tcur <= 4) %>% 
  group_by(year) %>% 
  filter(fit == max(fit)) %>% View()


daily_mbm_grid %>%
  unscale('bathy', ., resolution = 'fine') %>%
  filter(Species_code == 'CoMu') %>% 
  dplyr::select(zone, bathy) %>% 
  distinct()

# dth
CMFine2_dth %>%
  unscale('dth', ., resolution = 'fine') %>%
  filter(dth >= 2.3) %>%
  group_by(year) %>% 
  filter(fit == max(fit)) %>% View()


daily_mbm_grid %>%
  unscale('bathy', ., resolution = 'fine') %>%
  filter(Species_code == 'CoMu') %>% 
  dplyr::select(zone, bathy) %>% 
  distinct()

# variability explianed by null model for gulls
CoMu_year <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'CoMu') %>% 
  mgcv::gam(Count ~ s(year, bs="re"),
            data = .,
            offset = log(Effort_sqkm),
            family = 'nb')
summary(CoMu_year)


# Harbor Seal -------------------------------------------------------------

# dist
HSFine2_dist %>% 
  unscale('dist', ., resolution = 'fine') %>%
  filter(dist >= 500 & dist <= 900) %>%
  group_by(year) %>% 
  filter(fit == min(fit)) %>% View()

daily_mbm_grid %>%
  unscale('dist', ., resolution = 'fine') %>%
  filter(Species_code == 'HSeal') %>% 
  dplyr::select(zone, dist) %>% 
  distinct()

# sst
HSFine2_sst %>% 
  unscale('sst', ., resolution = 'fine') %>%
  # filter(dist >= 500 & dist <= 900) %>%
  group_by(year) %>% 
  filter(fit == max(fit)) %>% View()

# assess random effects
HSeal_year <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(year, bs="re"),
            data = .,
            family = 'binomial')

summary(HSeal_year)

HSeal_cruise <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HSeal') %>% 
  mgcv::gam(PresAbs ~ s(cruise.gen, bs="re"),
            data = .,
            family = 'binomial')

summary(HSeal_cruise)


# HPorp -------------------------------------------------------------------

HPorp_year <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(year, bs="re"),
            data = .,
            family = 'binomial')

summary(HPorp_year)

HPorp_cruise <-
  daily_mbm_grid %>% 
  filter(Species_code == 'HPorp') %>% 
  mgcv::gam(PresAbs ~ s(cruise.gen, bs = "re"),
            data = .,
            family = 'binomial')

summary(HPorp_cruise)

