
# setup -------------------------------------------------------------------

library(tidyverse)

source('code/setup.R')
source('code/functions.R')



# Env. Condition Summary --------------------------------------------------

# phyto
daily_mbm_grid %>% 
  unscale('phyto', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = phyto)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = 'Chlorophyll Concentration (µmol/L)') +
  theme_classic()

# sst
daily_mbm_grid %>% 
  unscale('sst', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = sst)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = 'Sea-Surface Temperature (ºC)') +
  theme_classic()

# sdSST
daily_mbm_grid %>% 
  unscale('temp_sd', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = temp_sd)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = 'St. Deviation of SST (ºC)') +
  theme_classic()
  
# salt
daily_mbm_grid %>% 
  unscale('salt', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = salt)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = 'Sea Surface Salinity (PSU)') +
  theme_classic()

# dth
daily_mbm_grid %>% 
  unscale('dth', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = dth)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = '∆ Tide Height (m)') +
  theme_classic()

# Table 2 -----------------------------------------------------------------

daily_mbm_grid %>% 
  unscale('salt', data = .) %>% 
  group_by(year) %>% 
  summarize(
    SSS.min = min(salt, na.rm = T) %>% round(digits = 0),
    SSS.mean = mean(salt, na.rm = T) %>% round(digits = 0),
    SSS.max = max(salt, na.rm = T) %>% round(digits = 0),
    spacer = NA,
    Cholo.min = min(phyto, na.rm = T) %>% round(digits = 1),
    Cholo.mean = mean(phyto, na.rm = T) %>% round(digits = 1),
    Cholo.max = max(phyto, na.rm = T) %>% round(digits = 1),
    spacer = NA,
    SST.min = min(sst, na.rm = T) %>% round(digits = 1),
    SST.mean = mean(sst, na.rm = T) %>% round(digits = 1),
    SST.max = max(sst, na.rm = T) %>% round(digits = 1),
    spacer = NA,
    sdSST.min = min(temp_sd, na.rm = T) %>% round(digits = 2),
    sdSST.mean = mean(temp_sd, na.rm = T) %>% round(digits = 2),
    sdSST.max = max(temp_sd, na.rm = T) %>% round(digits = 2)
  ) %>% View()

# marine mammal abundance and number of cruises
mbm_data %>% 
  filter(Date <= lubridate::mdy('10-31-2021')) %>% 
  mutate(
    year = lubridate::year(Date),
    group = if_else(Species_code %in% c('GL', 'CoMu'),
                    'Seabird',
                    'MarMam')) %>% 
  group_by(year, group) %>% 
  summarize(
    Cruises = unique(Date) %>% length(),
    `Total Area` = (sum(Effort_sqkm)/2) %>% round(digits = 0),
    min.Density = min(Density) %>% round(digits = 2),
    Avg.Density = mean(Density) %>% round(digits = 2),
    max.Density = max(Density) %>% round(digits = 2)) %>% View()


# Table 5 -----------------------------------------------------------------


# summarize cross validation metrics
# year
LYO_metrics %>%
  group_by(species) %>%
  summarize(
    Dev.Expl = mean(Dev.Expl, na.rm = T),
    AUC = mean(AUC, na.rm = T),
    TSS = mean(TSS, na.rm = T),
    `Obs:Pred` = mean(`Obs:Pred`, na.rm = T),
    `sd.Obs:Pred` = mean(`sd.Obs:Pred`, na.rm = T),
    spearman.rho = mean(spearman.rho, na.rm = T),
    spearman.p = mean(spearman.p, na.rm = T)) %>% View()

# zone

LZO_metrics %>%
  group_by(species) %>%
  summarize(
    Dev.Expl = mean(Dev.Expl, na.rm = T),
    AUC = mean(AUC, na.rm = T),
    TSS = mean(TSS, na.rm = T),
    `Obs:Pred` = mean(`Obs:Pred`, na.rm = T),
    `sd.Obs:Pred` = mean(`sd.Obs:Pred`, na.rm = T),
    spearman.rho = mean(spearman.rho, na.rm = T),
    spearman.p = mean(spearman.p, na.rm = T)) %>% View()