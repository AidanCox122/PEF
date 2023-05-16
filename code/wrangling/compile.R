
# setup -------------------------------------------------------------------

library(tidyverse) 

# liveocean data
phyto <-
  read_csv('data/clean/phyto.csv')

temp <- 
  read_csv('data/clean/temp.csv')

salt <-
  read_csv('data/clean/salt.csv')


# cruise details ------------------------------------------------------------

# determine what grid cells fall within transect bounds
transect_zones <- 
  read_csv("data/transect_area.csv") %>%
  dplyr::select(id, name) %>% 
  rename(
    grid_id = id,
    zone = name
  )

cruises_ocean_time_all <- data.frame(
  Date = c("2017-10-03", "2017-10-10", "2017-10-24", "2017-10-31", "2017-11-07", "2017-11-16", "2018-10-02", "2018-10-11", "2018-10-16", "2018-10-23", "2018-10-30", "2018-11-06", "2019-10-02", "2019-10-19", "2019-10-23", "2019-10-30", "2019-11-05", "2020-10-06", "2020-10-07", "2020-10-15", "2020-10-20", "2020-10-26", "2020-11-02", "2020-11-10", "2021-10-07", "2021-10-12", "2021-10-19", "2021-10-26"),
  ocean_time = c(276, 283, 297, 304, 311, 320, 640, 649, 654, 661, 668, 675, 1005, 1022, 1026, 1033, 1034, 1375, 1376, 1384, 1389, 1395, 1402, 1410, 1741, 1746, 1753, 1760),
  cruise = c(seq(1,28))) %>% 
  mutate(Date = lubridate::date(Date))

# get environmental variables for each cruise date
## phytoplankton
LO_phyto <- 
  left_join(phyto, cruises_ocean_time_all) %>% 
  # select data from cruise dates
  filter(!is.na(cruise)) %>% 
  # find the average value in each grid cell on that day (4 LiveOcean points per grid)
  group_by(Date, grid_id) %>% 
  summarise(
    gridPhyto = mean(phyto, na.rm = TRUE)) %>% 
  # assign zones to each grid cell
  left_join(transect_zones, by = 'grid_id', relationship = 'many-to-many') %>% 
  # filter out grid cells with no zones
  filter(!is.na(zone)) %>% 
  # take the average daily value in each zone
  group_by(Date, zone) %>%
  summarise(
    zonePhyto = mean(gridPhyto)) %>% 
  # change zone to factor
  mutate(zone = factor(zone))

## temperature
# LO_temp <- 
  left_join(temp, cruises_ocean_time_all) %>% 
  # select data from cruise dates
  filter(!is.na(cruise)) %>% 
  # find the average value in each grid cell on that day (4 LiveOcean points per grid)
  group_by(Date, grid_id) %>% 
  summarise(
    gridTemp = mean(temp, na.rm = TRUE),
    temp_sd = sd(temp, na.rm = T)) %>% 
  # assign zones to each grid cell
  left_join(transect_zones, by = 'grid_id', relationship = 'many-to-many') %>% 
  # filter out grid cells with no zones
  filter(!is.na(zone)) %>% 
  # take the average daily value in each zone
  group_by(Date, zone) %>%
  summarise(
    zoneTemp = mean(gridTemp),
    temp_sd = mean(temp_sd)) %>% 
  # change zone to factor
  mutate(zone = factor(zone))
  
## salinity
LO_salt <-  
  left_join(salt, cruises_ocean_time_all) %>% 
  # select data from cruise dates
  filter(!is.na(cruise)) %>% 
  # find the average value in each grid cell on that day (4 LiveOcean points per grid)
  group_by(Date, grid_id) %>% 
  summarise(
    gridSalt = mean(salt, na.rm = TRUE)) %>% 
  # assign zones to each grid cell
  left_join(transect_zones, by = 'grid_id', relationship = 'many-to-many') %>% 
  # filter out grid cells with no zones
  filter(!is.na(zone)) %>% 
  # take the average daily value in each zone
  group_by(Date, zone) %>%
  summarise(
    zoneSalt = mean(gridSalt)) %>% 
  # change zone to factor
  mutate(zone = factor(zone))

LO_phyto <- LO_phyto %>% 
LO_temp <- LO_temp %>% filter(!is.na(cruise))
LO_salt <- LO_salt %>% filter(!is.na(cruise))