
# setup -------------------------------------------------------------------

library(tidyverse)

source('code/setup.R')
source('code/functions.R')



# Figure 1b --------------------------------------------------

# phyto
phyto <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
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
sst <- 
  daily_mbm_grid %>%
  filter(Species_code == 'GL') %>% 
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
sdSST <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
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
salt <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL') %>% 
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
dth <- 
  daily_mbm_grid %>% 
  filter(Species_code == 'GL' & zone == 1) %>% 
  unscale('dth', ., resolution = 'fine') %>% 
  dplyr::select(Date, zone, phyto, sst, temp_sd, salt, dth) %>%
  # pivot_longer(cols = phyto:dth, names_to = 'variable', values_to = 'Value') %>% 
  ggplot(aes(x = 1, y = dth)) +
  geom_jitter(alpha = 0.25) +
  geom_boxplot(fill = 'grey', alpha = 0.5) +
  scale_x_continuous(breaks = NULL) +
  labs(x = '', y = '∆ Tide Height (m)') +
  theme_classic()

## save fig 1b ----

env_comp_plots <-
  list(
    phyto,
    sst,
    sdSST,
    salt,
    dth) %>% 
  set_names(
    c(
      'phyto',
      'sst',
      'sdSST',
      'salt',
      'dth'))

for (x in names(env_comp_plots)) {
  fname <-
    paste0('products/figure1/raw/', Sys.Date(), '_', (x), '.tiff')
  ggsave(fname,
         env_comp_plots[[x]],
         device = 'tiff',
         width = 4.48,
         height = 3.36,
         dpi = 500,
         units = 'in')
  print(paste('Done with', x))}


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

