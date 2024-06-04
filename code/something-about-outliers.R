daily_mbm_grid %>% 
  mutate(Year = lubridate::year(Date)) %>% 
  group_by(Year, zone) %>% 
  summarize(
    avg.dth = mean(dth))

mydata_outlier <- 
  mydata_bydatezone %>%
  mutate(Year = lubridate::year(Date)) %>% 
  group_by(Year, Species_code) %>%
  filter(Count > quantile(Count, 0.75) + (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))) | 
           Count < quantile(Count, 0.25) - (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))))

mydata_regular <-
  mydata_bydatezone %>%
  mutate(Year = lubridate::year(Date)) %>% 
  group_by(Year, Species_code) %>%
  filter(Count <= quantile(Count, 0.75) + (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))) & 
           Count >= quantile(Count, 0.25) - (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))))

mydata_bydatezone <- 
  # find the maximum non-outlier count for each species in each year
  mydata_regular %>% 
  filter(Year >= 2017) %>% 
  group_by(Year, Species_code) %>%
  filter(Count == max(Count)) %>% 
  dplyr::select(Species_code, max.Count= Count, Year) %>% 
  distinct() %>% 
  # transform outlier values to highest non-outlier
  right_join(mydata_outlier, by = c('Species_code', 'Year')) %>% 
  dplyr::transmute(Species_code, Count = max.Count, Year, Date, Zone, Effort_sqkm) %>% 
  rbind(mydata_regular)
  
  
  mutate(Count = case_when(
    Species_code == 'GL' ~ max(mydata_regular %>% filter(Species_code == 'GL') %>% pull(Count)),
    Species_code == 'CoMu' ~ max(mydata_regular %>% filter(Species_code == 'CoMu') %>% pull(Count)),
    Species_code == 'HSeal' ~ max(mydata_regular %>% filter(Species_code == 'HSeal') %>% pull(Count)),
    Species_code == 'HPorp' ~ max(mydata_regular %>% filter(Species_code == 'HPorp') %>% pull(Count)))) %>% 
  rbind(mydata_regular)
