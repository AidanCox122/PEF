
# setup -------------------------------------------------------------------

library(tidyverse)
library(ncdf4)


# extract variables --------------------------------------------------------------

# Open the netcdf files containing model data and read in the layers containing variables of interest:

nc_data <- nc_open("data/cox_surf_2017.01.01_2021.10.30.nc")

lon <- ncvar_get(nc_data, "lon_rho")
lat <- ncvar_get(nc_data, "lat_rho")
t <- ncvar_get(nc_data, "ocean_time")

phyto.array <- ncvar_get(nc_data, "phytoplankton")
temp.array <- ncvar_get(nc_data, "temp")
sal.array <- ncvar_get(nc_data, "salt")

dim(phyto.array) # check the dimensions of the array, expected: 72(x) x 75(y) x 1764 days

# find out what value is input for misssing data
fillvalue <- ncatt_get(nc_data, "phytoplankton", "_FillValue")
phyto.array[phyto.array == fillvalue$value] <- NA

fillvalue <- ncatt_get(nc_data, "temp", "_FillValue")
temp.array[temp.array == fillvalue$value] <- NA

fillvalue <- ncatt_get(nc_data, "salt", "_FillValue")
sal.array[sal.array == fillvalue$value] <- NA

# close the NETCDF file
nc_close(nc_data) 

# Determine spatial points used by the model and where they fall in the study grid:

extraction_points <- data.frame (
  longitude = c(lon),
  latitude = c(lat)
)

# read in file with LiveOcean points by study grid cell
LiveOcean_grid <- 
  read.csv("data/LiveOcean_pts(grid).csv") %>% 
  dplyr::select(c(point, latitude, longitude, grid_id)) %>%
  rename(
  lat = latitude,
  lon = longitude
)

# extract variables -------------------------------------------------------


## phytoplankton -----------------------------------------------------------
r_brick <- brick(phyto.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# correct the orrientation of the raster brick
r_brick <- flip(t(r_brick), direction='y')

phyto_2017 <- data.frame()
phyto_2018 <- data.frame()
phyto_2019 <- data.frame()
phyto_2020 <- data.frame()
phyto_2021 <- data.frame()

for(i in 1:nrow(extraction_points)) {
  point_series <- raster::extract(r_brick, extraction_points[i,], method='simple')
  data <- data.frame(
    point = c(rep(i, times = 1764)),
    lat = rep(extraction_points[i,2], times = 1764), 
    lon = rep(extraction_points[i,1], times = 1764), 
    ocean_time = c(seq(1,1764)),
    phyto = c(point_series)
  )
  data17 <- data %>% filter(ocean_time >= 274 & ocean_time <= 334)
  data18 <- data %>% filter(ocean_time >= 639 & ocean_time <= 699)
  data19 <- data %>% filter(ocean_time >= 1004 & ocean_time <= 1064)
  data20 <- data %>% filter(ocean_time >= 1370 & ocean_time <= 1430)
  data21 <- data %>% filter(ocean_time >= 1735 & ocean_time <= 1764)
  phyto_2017 <- rbind(phyto_2017, data17)
  phyto_2018 <- rbind(phyto_2018, data18)
  phyto_2019 <- rbind(phyto_2019, data19)
  phyto_2020 <- rbind(phyto_2020, data20)
  phyto_2021 <- rbind(phyto_2021, data21)
  rm(data, data17, data18, data21, data19, data20)
  print(i)
}

# cleanup and formatting
phyto_2017 <- phyto_2017 %>% filter(!is.na(phyto))
phyto_2018 <- phyto_2018 %>% filter(!is.na(phyto))
phyto_2019 <- phyto_2019 %>% filter(!is.na(phyto))
phyto_2020 <- phyto_2020 %>% filter(!is.na(phyto))
phyto_2021 <- phyto_2021 %>% filter(!is.na(phyto))

phyto_2017 <- inner_join(phyto_2017, LiveOcean_grid, by = "point")
phyto_2018 <- inner_join(phyto_2018, LiveOcean_grid, by = "point")
phyto_2019 <- inner_join(phyto_2019, LiveOcean_grid, by = "point")
phyto_2020 <- inner_join(phyto_2020, LiveOcean_grid, by = "point")
phyto_2021 <- inner_join(phyto_2021, LiveOcean_grid, by = "point")

phyto_2017 <- phyto_2017 %>% filter(!is.na(grid_id))
phyto_2018 <- phyto_2018 %>% filter(!is.na(grid_id))
phyto_2019 <- phyto_2019 %>% filter(!is.na(grid_id))
phyto_2020 <- phyto_2020 %>% filter(!is.na(grid_id))
phyto_2021 <- phyto_2021 %>% filter(!is.na(grid_id))

rm(r_brick) #cleanup


## sea-surface temperature ---------------------------------------------------------------------

temp_brick <- brick(temp.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# correct the orrientation of the raster brick
temp_brick <- flip(t(temp_brick), direction='y')

temp_2017 <- data.frame()
temp_2018 <- data.frame()
temp_2019 <- data.frame()
temp_2020 <- data.frame()
temp_2021 <- data.frame()

for(i in 1:nrow(extraction_points)) {
  point_series <- raster::extract(temp_brick, extraction_points[i,], method='simple')
  data <- data.frame(
    point = c(rep(i, times = 1764)),
    lat = rep(extraction_points[i,2], times = 1764), 
    lon = rep(extraction_points[i,1], times = 1764), 
    ocean_time = c(seq(1,1764)),
    temp = c(point_series)
  )
  data17 <- data %>% filter(ocean_time >= 274 & ocean_time <= 334)
  data18 <- data %>% filter(ocean_time >= 639 & ocean_time <= 699)
  data19 <- data %>% filter(ocean_time >= 1004 & ocean_time <= 1064)
  data20 <- data %>% filter(ocean_time >= 1370 & ocean_time <= 1430)
  data21 <- data %>% filter(ocean_time >= 1735 & ocean_time <= 1764)
  temp_2017 <- rbind(temp_2017, data17)
  temp_2018 <- rbind(temp_2018, data18)
  temp_2019 <- rbind(temp_2019, data19)
  temp_2020 <- rbind(temp_2020, data20)
  temp_2021 <- rbind(temp_2021, data21)
  rm(data, data17, data18, data21, data19, data20)
  print(i)
}

# cleanup and formatting
temp_2017 <- temp_2017 %>% filter(!is.na(temp))
temp_2018 <- temp_2018 %>% filter(!is.na(temp))
temp_2019 <- temp_2019 %>% filter(!is.na(temp))
temp_2020 <- temp_2020 %>% filter(!is.na(temp))
temp_2021 <- temp_2021 %>% filter(!is.na(temp))

temp_2017 <- inner_join(temp_2017, LiveOcean_grid, by = "point")
temp_2018 <- inner_join(temp_2018, LiveOcean_grid, by = "point")
temp_2019 <- inner_join(temp_2019, LiveOcean_grid, by = "point")
temp_2020 <- inner_join(temp_2020, LiveOcean_grid, by = "point")
temp_2021 <- inner_join(temp_2021, LiveOcean_grid, by = "point")

temp_2017 <- temp_2017 %>% filter(!is.na(grid_id))
temp_2018 <- temp_2018 %>% filter(!is.na(grid_id))
temp_2019 <- temp_2019 %>% filter(!is.na(grid_id))
temp_2020 <- temp_2020 %>% filter(!is.na(grid_id))
temp_2021 <- temp_2021 %>% filter(!is.na(grid_id))

rm(temp_brick)


## salinity ----------------------------------------------------------------

## Data Extraction - SALINITY ##
sal_brick <- brick(sal.array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# correct the orrientation of the raster brick
sal_brick <- flip(t(sal_brick), direction='y')

sal_2017 <- data.frame()
sal_2018 <- data.frame()
sal_2019 <- data.frame()
sal_2020 <- data.frame()
sal_2021 <- data.frame()

for(i in 1:nrow(extraction_points)) {
  point_series <- raster::extract(sal_brick, extraction_points[i,], method='simple')
  data <- data.frame(
    point = c(rep(i, times = 1764)),
    lat = rep(extraction_points[i,2], times = 1764), 
    lon = rep(extraction_points[i,1], times = 1764), 
    ocean_time = c(seq(1,1764)),
    salt = c(point_series)
  )
  data17 <- data %>% filter(ocean_time >= 274 & ocean_time <= 334)
  data18 <- data %>% filter(ocean_time >= 639 & ocean_time <= 699)
  data19 <- data %>% filter(ocean_time >= 1004 & ocean_time <= 1064)
  data20 <- data %>% filter(ocean_time >= 1370 & ocean_time <= 1430)
  data21 <- data %>% filter(ocean_time >= 1735 & ocean_time <= 1764)
  sal_2017 <- rbind(sal_2017, data17)
  sal_2018 <- rbind(sal_2018, data18)
  sal_2019 <- rbind(sal_2019, data19)
  sal_2020 <- rbind(sal_2020, data20)
  sal_2021 <- rbind(sal_2021, data21)
  rm(data, data17, data18, data21, data19, data20)
  print(i)
}

# cleanup and formatting
sal_2017 <- sal_2017 %>% filter(!is.na(salt))
sal_2018 <- sal_2018 %>% filter(!is.na(salt))
sal_2019 <- sal_2019 %>% filter(!is.na(salt))
sal_2020 <- sal_2020 %>% filter(!is.na(salt))
sal_2021 <- sal_2021 %>% filter(!is.na(salt))

sal_2017 <- inner_join(sal_2017, LiveOcean_grid, by = "point")
sal_2018 <- inner_join(sal_2018, LiveOcean_grid, by = "point")
sal_2019 <- inner_join(sal_2019, LiveOcean_grid, by = "point")
sal_2020 <- inner_join(sal_2020, LiveOcean_grid, by = "point")
sal_2021 <- inner_join(sal_2021, LiveOcean_grid, by = "point")

sal_2017 <- sal_2017 %>% filter(!is.na(grid_id))
sal_2018 <- sal_2018 %>% filter(!is.na(grid_id))
sal_2019 <- sal_2019 %>% filter(!is.na(grid_id))
sal_2020 <- sal_2020 %>% filter(!is.na(grid_id))
sal_2021 <- sal_2021 %>% filter(!is.na(grid_id))

rm(sal_brick)


# compile
phyto <- rbind(phyto_2017, phyto_2018, phyto_2019, phyto_2020, phyto_2021)
temp <- rbind(temp_2017, temp_2018, temp_2019, temp_2020, temp_2021)
salt <- rbind(sal_2017, sal_2018, sal_2019, sal_2020, sal_2021)

write_csv(phyto, 'data/clean/phyto.csv')
write_csv(temp, 'data/clean/temp.csv')
write_csv(salt, 'data/clean/salt.csv')

