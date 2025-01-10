
# HSeal -------------------------------------------------------------------
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

# HPorp -------------------------------------------------------------------
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


# GL ----------------------------------------------------------------------

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
  mgcv::gam(Count ~ s(tcur, bs = 'bs', m=c(3,1), k = 5) + s(dth, bs = 'bs', m=c(3,1)) + s(sst, bs = 'bs', m=c(3,1)) + s(year, bs="re"),
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

# CoMu --------------------------------------------------------------------

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



