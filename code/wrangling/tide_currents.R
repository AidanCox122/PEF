# DOWNLOADING TIDAL DATA 


# setup -------------------------------------------------------------------

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(sf)

world <- ne_countries(scale = "medium", returnclass = "sf")

# Station info -------------------------------------------------------------------

# Here is a list of the tide height and current stations within my study area and their locations:
  
## Current Prediction Stations for Consideration
current_stations <- data.frame(
  station = c("PCT1966", "PCT1971", "PUG1728", "PUG1727", "PUG1702", "PUG1730", "PCT2006", "PUG1729", "PCT2026", "PUG1731", "PUG1733", "PCT2046", "PUG1732", "PUG1704", "PUG1705", "PUG1706", "PCT2071", "PCT2076", "PCT2126", "PUG1707", "PUG1708", "PUG1712", "PCT1416", "PUG1742", "PUG1703", "PCT2191", "PUG1746", "PUG1745", "PUG1723", "PUG1721", "PUG1720", "PUG1719", "PUG1715", "PUG1722", "PUG1744", "PUG1724", "PUG1718", "PCT2266", "PUG1716", "PCT2281"),
  location = c("Iceberg Pass", "Colville Island", "Point Colville", "Lawson Reef", "Rosario Strait", "Lopez Pass", "Burrows Bay", "Belle Rock", "Green Point", "Fontleroy Light", "Thatcher Pass", "Frost Willow Island", "Strawberry Island", "Peavine Pass", "Obstruction Pass", "Peapod Rocks", "Barnes Island", "Raccoon Island", "Towhead Island", "Sinclair Island", "Lawrence Point", "Parker Reef", "Cattle Point", "Cattle Point 2", "SJC South", "King's Point", "Pear Point", "Point George", "Upright Channel", "Wasp Passage", "Spring Passage", "Spieden Channel", "President Channel", "Harney Passage", "Discovery Island", "Lime Kiln", "Kellett Bluff", "John's Island", "Waldron Island", "Point Hammond" ),
  latitude = c(48.3833, 48.4000, 48.4181, 48.4125, 48.4581, 48.4797, 48.4628, 48.4968, 48.5047, 48.5216, 48.5274, 48.5392, 48.5610, 48.5871, 48.6033, 48.6224, 48.6858, 48.6122, 48.6442, 48.6794, 48.7326, 48.4000, 48.3840, 48.4344, 48.4610, 48.4833, 48.5114, 48.5567, 48.5538, 48.5925, 48.6115, 48.6278, 48.6734, 48.5897, 48.4521, 48.4980, 48.5887, 48.6833, 48.7042, 48.7320 ),
  longitude = c(-122.9167, -122.8167, -122.7812, -122.7403, -122.7501, - 122.8189, -122.6828, -122.7308, -122.7062, -122.7707, -122.8040, -122.8308, -122.7543, -122.8193, -122.8127, -122.7476, -122.7888, -122.7022, -122.6587, -122.7147, -122.8864, -123.0000, -123.0157, -122.9466, -122.9520, -122.9558, -122.9529, -122.9985, -122.9226, -122.9896, -123.0341, -123.1116, -123.0060, -122.9217, -123.1554, -123.1599, -123.2258, -123.1500, -123.1048, -123.0253)
)

## Map Station Locations
ggplot() +
  geom_sf(data = world, color = "black") +
  coord_sf(xlim = c(-124, -122), ylim = c(48,49), expand = FALSE) +
  
  geom_point(data = current_stations, aes(longitude, latitude, color = location)) +
  
  xlab("Longitude") +
  ylab("Latitude") +
  theme_minimal() +
  theme(legend.position = "none")


# TIDE CURRENT ------------------------------------------------------------
## 2019 ##
fList <- list.files("data/Tides/2019/", full.names = T) # update with file path from chunk 4

current_amplitudes_2019 <- data.frame() 

for(i in 1:length(fList)) {
  currents <- 
    read_csv(fList[i]) %>%
    # split DateTime into untidy but useful subdivisions
    mutate(
    DateTime = `Date_Time (LST/LDT)`) %>% 
    separate(DateTime, into = c("Date", "Time"), sep = "([ ])") %>%
    separate(Date, into = c("Year", "Month", "Day"), sep = "([-])")%>% 
    # select only data from October and November
    filter(Month == 10 | Month == 11) %>%
    # find the maximum flow and ebb speeds for each date
    filter(`Speed (knots)` != "-") %>% 
    mutate(`Speed (knots)` = as.numeric(`Speed (knots)`)) %>%
    group_by(Month, Day, Event) %>%
    top_n(1, abs(`Speed (knots)`)) %>%
    # select important variables
    dplyr::select(
    Year, Month, Day, Event, `Speed (knots)`) %>% 
    distinct() %>% 
    # calculate the amplitude of tidal current change from max flood to max ebb
    mutate(speed = abs(`Speed (knots)`)) %>%
    group_by(Month, Day) %>%
    mutate(Amplitude = speed + lead(speed),
           Station = fList[i]) %>% 
    ungroup() %>% 
    separate(Station,
             into = c("Trash1", "Trash2", "Trash3", "Trash4", "Trash5", "Trash6", "Trash7", "Trash8", "Station"),
             sep = "([/])") %>% 
    separate(Station,
             into = c("Station", "Trash9", "Trash10", "Trash11"),
             sep = "([_])")%>%
    filter(Event == 'flood' & Amplitude != 0) %>% 
    dplyr::select(-c(`Speed (knots)`)) 
  name <- paste(currents[1,1], currents[1,'Station'], sep = "_")
  #name <- paste(currents[1,1], currents[1,5], sep = "_")
  #write_csv(currents, paste("~/Desktop/PEF/2019/", name, ".csv", sep = ""))
  current_amplitudes_2019 <- rbind(currents, current_amplitudes_2019) # update with year x2
  rm(currents, name)
  print(i)
}
current_amplitudes_2019$Year <- "2019" 

## 2020 ##
fList <- list.files("~/Desktop/PEF/Tides/2020/", full.names = T) # update with file path from chunk 4

current_amplitudes_2020 <- data.frame() 

for(i in 1:length(fList)) {
  currents <- read_csv(fList[i]) 
  currents <- currents %>% mutate(
    DateTime = `Date_Time (LST/LDT)`
  )
  currents <- currents %>% separate(DateTime, into = c("Date", "Time"), sep = "([ ])")
  currents <- currents %>% separate(Date, into = c("Year", "Month", "Day"), sep = "([-])")
  currents <- currents %>% filter(Month == 10 | Month == 11) # select only data from October and November
  currents <- currents %>% filter(`Speed (knots)` != "-")
  currents$`Speed (knots)` <- as.numeric(currents$`Speed (knots)`)
  currents <- currents %>% group_by(Month, Day, Event) %>%
    top_n(1, abs(`Speed (knots)`))
  currents <- currents %>% dplyr::select(
    Year, Month, Day, Event, `Speed (knots)`
  )
  currents <- distinct(currents)
  currents$`Speed (knots)` <- abs(currents$`Speed (knots)`)
  currents <- currents %>% group_by(Month, Day) %>%
    mutate(Amplitude = `Speed (knots)` + lag(`Speed (knots)`))
  currents$Station <- fList[i]
  currents <- currents %>% separate(Station, into = c("Trash1", "Trash2", "Trash3", "Trash4", "Trash5", "Trash6", "Trash7", "Trash8", "Station"), sep = "([/])")
  currents <- currents %>% separate(Station, into = c("Station", "Trash9", "Trash10", "Trash11"), sep = "([_])")
  currents <- currents %>% dplyr::select(Year, Month, Day, Amplitude, Station)
  currents <- currents %>% filter(Amplitude != 0)
  name <- paste(currents[1,1], currents[1,5], sep = "_")
  #name <- paste(currents[1,1], currents[1,5], sep = "_")
  #write_csv(currents, paste("~/Desktop/PEF/2020/", name, ".csv", sep = ""))
  currents <- ungroup(currents)
  current_amplitudes_2020 <- rbind(currents, current_amplitudes_2020) # update with year x2
  rm(currents, name)
  print(i)
}
current_amplitudes_2020$Year <- "2020" 

## 2021 ##
fList <- list.files("~/Desktop/PEF/Tides/2021/", full.names = T) # update with file path from chunk 4

current_amplitudes_2021 <- data.frame() 

for(i in 1:length(fList)) {
  currents <- read_csv(fList[i]) 
  currents <- currents %>% mutate(
    DateTime = `Date_Time (LST/LDT)`
  )
  currents <- currents %>% separate(DateTime, into = c("Date", "Time"), sep = "([ ])")
  currents <- currents %>% separate(Date, into = c("Year", "Month", "Day"), sep = "([-])")
  currents <- currents %>% filter(Month == 10 | Month == 11) # select only data from October and November
  currents <- currents %>% filter(`Speed (knots)` != "-")
  currents$`Speed (knots)` <- as.numeric(currents$`Speed (knots)`)
  currents <- currents %>% group_by(Month, Day, Event) %>%
    top_n(1, abs(`Speed (knots)`))
  currents <- currents %>% dplyr::select(
    Year, Month, Day, Event, `Speed (knots)`
  )
  currents <- distinct(currents)
  currents$`Speed (knots)` <- abs(currents$`Speed (knots)`)
  currents <- currents %>% group_by(Month, Day) %>%
    mutate(Amplitude = `Speed (knots)` + lag(`Speed (knots)`))
  currents$Station <- fList[i]
  currents <- currents %>% separate(Station, into = c("Trash1", "Trash2", "Trash3", "Trash4", "Trash5", "Trash6", "Trash7", "Trash8", "Station"), sep = "([/])")
  currents <- currents %>% separate(Station, into = c("Station", "Trash9", "Trash10", "Trash11"), sep = "([_])")
  currents <- currents %>% dplyr::select(Year, Month, Day, Amplitude, Station)
  currents <- currents %>% filter(Amplitude != 0)
  name <- paste(currents[1,1], currents[1,5], sep = "_")
  #name <- paste(currents[1,1], currents[1,5], sep = "_")
  #write_csv(currents, paste("~/Desktop/PEF/2021/", name, ".csv", sep = ""))
  currents <- ungroup(currents)
  current_amplitudes_2021 <- rbind(currents, current_amplitudes_2021) # update with year x2
  rm(currents, name)
  print(i)
}

current_amplitudes_2021$Year <- "2021" 


# TIDE HEIGHT -------------------------------------------------------------

# 2017
delta_2017 <-
  # read tide height data from november
  read_csv("data/Tide_Height/2017/2017_9449880_novacdc28a1638b.csv") %>% 
  # add data from october
  rbind(
    read_csv("data/Tide_Height/2017/2017_9449880_octacdc6c288d25.csv")) %>%
  # duplicate Date Time column and format correctly
  mutate(
    dup = `Date Time`) %>% 
  separate(dup, into = c("Date", "Time"), sep = "([ ])") %>% 
  mutate(
    Date = lubridate::ymd(Date))

# 2018
delta_2018 <-
  # read tide height predictions from november
  read_csv("data/Tide_Height/2018/2018_9449880_octacdc74ea0e9.csv") %>% 
  # add predictions from october
  rbind(
    read_csv("data/Tide_Height/2018/2018_9449880_novacdc5b6a46ec.csv")) %>% 
  # duplicate and format date 
  mutate(
    dup = `Date Time`) %>% 
    separate(dup, into = c("Date", "Time"), sep = "([ ])") %>% 
    mutate(
      Date = lubridate::ymd(Date))

# 2019
delta_2019 <-
  # read tide height predictions from november
  read_csv("data/Tide_Height/2019/2019_9449880_octacdca1cd424.csv") %>% 
  # add predictions from october
  rbind(
    read_csv("data/Tide_Height/2019/2019_9449880_novacdc5698ead8.csv")) %>% 
  # duplicate and format date 
  mutate(
    dup = `Date Time`) %>% 
  separate(dup, into = c("Date", "Time"), sep = "([ ])") %>% 
  mutate(
    Date = lubridate::ymd(Date))

# 2020
delta_2020 <-
  # read tide height predictions from november
  read_csv("data/Tide_Height/2020/2020_9449880_octacdcebe78ce.csv") %>% 
  # add predictions from october
  rbind(
    read_csv("data/Tide_Height/2020/2020_9449880_novacdc1a32ace5.csv")) %>% 
  # duplicate and format date 
  mutate(
    dup = `Date Time`) %>% 
  separate(dup, into = c("Date", "Time"), sep = "([ ])") %>% 
  mutate(
    Date = lubridate::ymd(Date))

# 2021
delta_2021 <-
  # read tide height predictions from november
  read_csv("data/Tide_Height/2021/2021_9449880_octc51f32ede5b0.csv") %>% 
  # add predictions from october
  rbind(
    read_csv("data/Tide_Height/2021/2021_9449880_novc51f77960932.csv")) %>% 
  # duplicate and format date 
  mutate(
    dup = `Date Time`) %>% 
  separate(dup, into = c("Date", "Time"), sep = "([ ])") %>% 
  mutate(
    Date = lubridate::ymd(Date))


delta.tide.height <- rbind(delta_2017, delta_2018, delta_2019, delta_2020, delta_2021)
rm(delta_2017, delta_2018, delta_2019, delta_2020, delta_2021)

# find lowest low and highest high for each day
delta.tide.height <- 
  delta.tide.height %>%
  group_by(Date) %>% 
  summarize(HH = max(Prediction),
            LL = min(Prediction)) %>% 
  mutate(dth = (HH - LL))




