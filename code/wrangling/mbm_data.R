# FORMAT MARINE BIRD AND MAMMAL DATA

# setup -------------------------------------------------------------------

library(tidyverse) 
library(readxl)


# formatting --------------------------------------------------------------

# read information on species codes, names and groups
species_names <- read_excel("data/MBM_data/MBM_Data2008-2021_MASTER.xlsx", range = "metadata!A10:f45")

# read data
byrecord <- read_excel("data/MBM_data/MBM_Data2008-2021_MASTER.xlsx", range = "countbyrecord!A1:AN30000")

# select rows with data
byrecord_clean <-
  byrecord %>% 
  filter(!is.na(Date) & Date>0) %>%    # remove blank lines: n = 18,488
  # standardize time and date and year
  mutate(Time = hms::as_hms(Time),
         Date = ymd(Date),
         Year = year(Date))


### Sampling effort first ###
# read sampling effort (number of times a zone is sampled during a cruise) from a file
effort <- read_excel("data/MBM_data/MBM_Data2008-2021_MASTER.xlsx", range = "sampling effort!A14:g150")
effort$Date <- ymd(effort$Date)  # date format
# select rows with data
effort %>% 
  filter(Date>0) ->  # remove blank lines
  effort             # n = 95 records

# read area sizes (square kilometers) for each zone
area_size <- read_excel("data/MBM_data/MBM_Data2008-2021_MASTER.xlsx", range = "sampling effort!A2:b8")

# create tibble to translate zone naming, e.g., Zone_1 to 1
zone_names <- 
  tibble(
    Zone_name = c("Zone_1","Zone_2","Zone_3","Zone_4","Zone_5","Zone_6"), 
    Zone = seq(1,6)
  )

# reformat effort 
effort %>%
  gather(key=Zone_name,value=Count,"Zone_1":"Zone_6") %>%
  left_join(zone_names,by = "Zone_name") %>%
  left_join(area_size,by = "Zone") ->  
  effort_long                                 # n = 534 records 

# compute area sampled for each date and zone
effort_long$Effort_sqkm <- effort_long$Count * effort_long$Size_sqkm

### Number of MBM is next ###
# sum species counts by Date and Zone 
byrecord_clean %>%
  dplyr::select(-c(Year)) %>%                   # drop year from summing
  group_by(Date,Zone) %>%                       # group by date and zone
  summarize_if(is.numeric, sum, na.rm=TRUE) ->  # sum by date
  byrecord_sumbyzone                              # n = 568 records 

# join effort data and species count data by Date and Zone
effort_long %>%
  left_join(byrecord_sumbyzone, by = c("Date","Zone")) %>%
  dplyr::select(c("Date","Zone","Effort_sqkm","GL","CoMu", 'HSeal', 'HPorp')) %>%
  gather(key=Species_code,value=Count,"GL":"HPorp") ->
  mydata_bydatezone   # n = 2270

# convert NA Count to 0 (these are sampled zones with no MBM observed)
mydata_bydatezone$Count[is.na(mydata_bydatezone$Count)] <- 0

# remove unsampled zones (Effort = 0)
mydata_bydatezone %>%
  dplyr::filter(Effort_sqkm > 0) ->
  mydata_bydatezone

# transform high outliers into highest regular count
mydata_outlier <- 
  mydata_bydatezone %>%
  group_by(Species_code) %>%
  filter(Count > quantile(Count, 0.75) + (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))) | 
           Count < quantile(Count, 0.25) - (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))))

mydata_regular <-
  mydata_bydatezone %>%
  group_by(Species_code) %>%
  filter(Count <= quantile(Count, 0.75) + (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))) & 
           Count >= quantile(Count, 0.25) - (1.5 * (quantile(Count,0.75) - quantile(Count,0.25))))

mydata_bydatezone <- 
  mydata_outlier %>% 
  group_by(Species_code) %>% 
  mutate(Count = case_when(
    Species_code == 'GL' ~ max(mydata_regular %>% filter(Species_code == 'GL') %>% pull(Count)),
    Species_code == 'CoMu' ~ max(mydata_regular %>% filter(Species_code == 'CoMu') %>% pull(Count)),
    Species_code == 'HSeal' ~ max(mydata_regular %>% filter(Species_code == 'HSeal') %>% pull(Count)),
    Species_code == 'HPorp' ~ max(mydata_regular %>% filter(Species_code == 'HPorp') %>% pull(Count)))) %>% 
  rbind(mydata_regular)

# compute species density by Date and Zone
mydata_bydatezone$Density <- round(mydata_bydatezone$Count / mydata_bydatezone$Effort_sqkm,2)

write_csv(mydata_bydatezone, 'data/clean/mbm_master.csv')
