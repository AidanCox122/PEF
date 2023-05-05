# DOWNLOADING TIDAL DATA 


# Station info -------------------------------------------------------------------

## EXTRACTING TIDAL CURRENT DATA

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


# extract raw data --------------------------------------------------------




