
# setup -------------------------------------------------------------------

library(tidyverse)
library(curl)

current_stations <- data.frame(
  station = c("PCT1966", "PCT1971", "PUG1728", "PUG1727", "PUG1702", "PUG1730", "PCT2006", "PUG1729", "PCT2026", "PUG1731", "PUG1733", "PCT2046", "PUG1732", "PUG1704", "PUG1705", "PUG1706", "PCT2071", "PCT2076", "PCT2126", "PUG1707", "PUG1708", "PUG1712", "PCT1416", "PUG1742", "PUG1703", "PCT2191", "PUG1746", "PUG1745", "PUG1723", "PUG1721", "PUG1720", "PUG1719", "PUG1715", "PUG1722", "PUG1744", "PUG1724", "PUG1718", "PCT2266", "PUG1716", "PCT2281"),
  location = c("Iceberg Pass", "Colville Island", "Point Colville", "Lawson Reef", "Rosario Strait", "Lopez Pass", "Burrows Bay", "Belle Rock", "Green Point", "Fontleroy Light", "Thatcher Pass", "Frost Willow Island", "Strawberry Island", "Peavine Pass", "Obstruction Pass", "Peapod Rocks", "Barnes Island", "Raccoon Island", "Towhead Island", "Sinclair Island", "Lawrence Point", "Parker Reef", "Cattle Point", "Cattle Point 2", "SJC South", "King's Point", "Pear Point", "Point George", "Upright Channel", "Wasp Passage", "Spring Passage", "Spieden Channel", "President Channel", "Harney Passage", "Discovery Island", "Lime Kiln", "Kellett Bluff", "John's Island", "Waldron Island", "Point Hammond" ),
  latitude = c(48.3833, 48.4000, 48.4181, 48.4125, 48.4581, 48.4797, 48.4628, 48.4968, 48.5047, 48.5216, 48.5274, 48.5392, 48.5610, 48.5871, 48.6033, 48.6224, 48.6858, 48.6122, 48.6442, 48.6794, 48.7326, 48.4000, 48.3840, 48.4344, 48.4610, 48.4833, 48.5114, 48.5567, 48.5538, 48.5925, 48.6115, 48.6278, 48.6734, 48.5897, 48.4521, 48.4980, 48.5887, 48.6833, 48.7042, 48.7320 ),
  longitude = c(-122.9167, -122.8167, -122.7812, -122.7403, -122.7501, - 122.8189, -122.6828, -122.7308, -122.7062, -122.7707, -122.8040, -122.8308, -122.7543, -122.8193, -122.8127, -122.7476, -122.7888, -122.7022, -122.6587, -122.7147, -122.8864, -123.0000, -123.0157, -122.9466, -122.9520, -122.9558, -122.9529, -122.9985, -122.9226, -122.9896, -123.0341, -123.1116, -123.0060, -122.9217, -123.1554, -123.1599, -123.2258, -123.1500, -123.1048, -123.0253)
)

# download data -----------------------------------------------------------

#The parts of the url to be assembled:
url1 = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date="
url2 = "&end_date="
url3 = "&station=" #return stationId
url4 = "&product=currents" #return product
url5 = "&datum=mllw" #return datum
url6 = "&units=metric" #return units
url7 = "&time_zone=lst" #return time zone
url8 = "&application=univer_washington" #return application
url9 = "&format=csv" #return format
### 

stations = c("PCT1966", "PCT1971", "PUG1728", "PUG1727", "PUG1702", "PUG1730", "PCT2006", "PUG1729", "PCT2026", "PUG1731", "PUG1733", "PCT2046", "PUG1732", "PUG1704", "PUG1705", "PUG1706", "PCT2071", "PCT2076", "PCT2126", "PUG1707", "PUG1708", "PUG1712", "PCT1416", "PUG1742", "PUG1703", "PCT2191", "PUG1746", "PUG1745", "PUG1723", "PUG1721", "PUG1720", "PUG1719", "PUG1715", "PUG1722", "PUG1744", "PUG1724", "PUG1718", "PCT2266", "PUG1716", "PCT2281")

## 2019
dir <- "data/Tides/2019" # set file-path

# 2019 October - November
for(i in 1:length(stations)) {
  begin_date <- 20191020
  end_date <- 20191120
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2019", stations[i], "oct-nov", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(url = urltotal, destfile = tmp)
}

# 2019 November - December
for(i in 1:length(stations)) {
  begin_date <- 20191121
  end_date <- 20191210
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2019", stations[i], "nov-dec", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(urltotal, tmp)
}

## 2020
dir <- "data/raw/Tides/2020"

# 2020 October-November
for(i in 1:length(stations)) {
  begin_date <- 20201001
  end_date <- 20201101
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2020", stations[i], "oct-nov", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(urltotal, tmp)
}

# 2020 November-December
for(i in 1:length(stations)) {
  begin_date <- 20201101
  end_date <- 20201201
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2020", stations[i], "nov-dec", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(urltotal, tmp)
}

## 2021
dir <- "data/raw/Tides/2021"

# 2021 October-November
for(i in 1:length(stations)) {
  begin_date <- 20211001
  end_date <- 20211101
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2021", stations[i], "oct-nov", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(urltotal, tmp)
}

# 2021 November-December
for(i in 1:length(stations)) {
  begin_date <- 20211101
  end_date <- 20211201
  station <- stations[i]
  urltotal <- paste(url1,begin_date,url2,end_date,url3,station,url4,url6,url7,url8,url9,sep ="")
  filename <- paste("2021", stations[i], "nov-dec", sep = "_")
  tmp <- tempfile(pattern = filename, tmpdir = dir, fileext = ".csv")
  curl_download(urltotal, tmp)
}
