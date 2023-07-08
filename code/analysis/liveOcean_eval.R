## THIS CODE CHUNK EVALUATES THE ACCURACY OF ENVIRONMENTAL DATA FROM 
## LIVEOCEAN BY COMPARING IT TO MEASURED ENV. DATA FROM PEF CTD CASTS


# setup -------------------------------------------------------------------

library(tidyverse)

# read in data from PEF CTD casts
# note: this data was transcribed by hand from raw record
pef_ctd <- 
  read_csv('data/PEF_CTD.csv')

# read LiveOcean data
temp <- 
  read_csv('data/clean/temp.csv') %>% 
  # convert lat.x and lon.x to characters for join key
  mutate(lat.x = as.character(lat.x),
         lon.x = as.character(lon.x))

salt <- 
  read_csv('data/clean/salt.csv') %>% 
  # convert lat.x and lon.x to characters for join key
  mutate(lat.x = as.character(lat.x),
         lon.x = as.character(lon.x))

# location of 4 LiveOcean points nearest ctd stations
model_pt <- data.frame(
  Station = c("N", "N", "N", "N", "S", "S", "S", "S"),
  lat.y = c(48.5825895797052,48.5825895797052,48.5870894726895,48.5870894726895, 48.4205934322702, 48.4205934322702,48.4160935392859,48.4250933252545),
  lon.y = c(-123.038688775285,-123.045421593992,-123.038688775285,-123.045421593992,-122.944429313385,-122.937696494678,-122.944429313385,-122.944429313385)) %>% 
  mutate(
    lat.x = as.character(lat.y),
    lon.x = as.character(lon.y))

# format modelled data --------------------------------------------------------

# estimate modeled temperature at each CTD station
modTemp <- 
  # inner join by station points, should be 2192 obs. (274 days x 8 pts.)
  inner_join(model_pt, temp, by = c('lat.x', 'lon.x')) %>% 
  # create a data variable from oceantime
  mutate(Date = as.Date(ocean_time, origin = "2016-12-31")) %>% 
  # summarize temperature for 4 points around each station on each day
  group_by(Station, Date) %>% 
  summarise(ModTemp = mean(temp), .groups = 'drop')

# estimate modeled salinity at each CTD station
modSalt <- 
  # inner join by station points, should be 2192 obs. (274 days x 8 pts.)
  inner_join(model_pt, salt, by = c('lat.x', 'lon.x')) %>% 
  # create a data variable from oceantime
  mutate(Date = as.Date(ocean_time, origin = "2016-12-31")) %>% 
  # summarize temperature for 4 points around each station on each day
  group_by(Station, Date) %>% 
  summarise(ModSalt = mean(salt), .groups = 'drop')


# combine modeled data with observed data
LivOce_eval <-
  inner_join(pef_ctd, modTemp) %>% 
  inner_join(modSalt)

# perform a paired t-test
# TEMPERATURE #
temp_test <- 
  t.test(
    LivOce_eval$ObsTemp,
    LivOce_eval$ModTemp,
    paired = TRUE,
    alternative = "two.sided")

print(temp_test)

# SALINITY #
salt_test <- 
  t.test(
    LivOce_eval$ObsSalt,
    LivOce_eval$ModSalt,
    paired = TRUE,
    alternative = "two.sided")

print(salt_test)



