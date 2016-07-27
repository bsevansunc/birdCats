# =================================================================================*
# ---- SET-UP ----
# =================================================================================*

library(tigris)
library(acs)
library(stringr) # to pad fips codes
library(leaflet)
library(rgdal)
library(dplyr)
library(sp)

options(stringsAsFactors = FALSE)

# =================================================================================*
# ---- GET CENSUS DATA ----
# =================================================================================*

# Add api key for census data:

api.key.install(key = 'ec7729865a6a6bdf5d1ce20999a4788da99b6c5e')

# Get census tracts for each state:

tractsMD <- tracts(state = 'MD')
tractsVA <- tracts(state = 'VA')
tractsDC <- tracts(state = 'DC')

# Bind census tracts into a single file:

tractsAll <- rbind_tigris(tractsMD, tractsVA, tractsDC)

# Take a look:

leaflet(tractsAll) %>%
  addTiles() %>%
  addPolygons(popup = ~NAME)

# Make into a spatial file of each county and tract:

geo <- geo.make(state = c('DC', 'VA', 'MD'), county = '*', tract = '*')

# Fetch census data for tracts:
# This example relates to median household income. For more variables see: http://www2.census.gov/geo/tiger/TIGER_DP/2014ACS/Metadata/BG_METADATA_2014.txt

income <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                  table.number = "B19001", col.names = "pretty")

# Make into a data frame:

income_df <- data.frame(
  GEOID = paste0(str_pad(income@geography$state, 2, "left", pad="0"),
                 str_pad(income@geography$county, 3, "left", pad="0"),
                 str_pad(income@geography$tract, 6, "left", pad="0")),
  medianIncome = as.numeric(income1@estimate)
  )

# Merge census data with geographic dataset:

income_merged<- geo_join(tractsAll, income_df, "GEOID", "GEOID")

# ---------------------------------------------------------------------------------*
# ---- Plot census data ----
# ---------------------------------------------------------------------------------*

# Popup:

popup <- paste0("GEOID: ", income_merged$GEOID, "<br>", "Median household income: ")

# Color palette:

pal <- colorNumeric(
  palette = "YlGnBu",
  domain = income_merged$medianIncome
)

map3<-leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = income_merged, 
              fillColor = ~pal(medianIncome), 
              # fillColor = ~pal(percent), 
              color = "#b2aeae", # you need to use hex colors
              fillOpacity = 0.7, 
              weight = 1, 
              smoothFactor = 0.2,
              popup = popup) %>%
  addLegend(pal = pal, 
            values = income_merged$medianIncome, 
            # values = income_merged$percent, 
            position = "bottomright", 
            title = "Median household<br>income ($)",
            # title = "Percent of Households<br>above $200k",
            labFormat = labelFormat(prefix = "")) 
map3

# =================================================================================*
# ---- CENSUS DATA AT CAT SAMPLING POINTS ----
# =================================================================================*

# Get nestwatch sampling coordinates (spatial coordinates are in UTM):

lc100 <- read.csv('lc100.csv')

# Project nestwatch data:

utm <- SpatialPoints(cbind(lc100$X,lc100$Y), proj4string=CRS("+proj=utm +zone=18"))

# Get a vector of cat sites:

catSites <- read.csv('catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT)) %>%
  select(SITE) %>%
  distinct %>%
  arrange(SITE) %>%
  .$SITE

# Transform sampling coordinates from UTM to lonlat and filter to cat sampling points:

lonlat <- spTransform(utm, CRS('+proj=longlat')) %>%
  data.frame %>%
  dplyr::rename(lon = coords.x1, lat = coords.x2) %>%
  cbind(lc100) %>%
  select(site, can, imp, lon, lat) %>%
  filter(site %in% catSites)

# Set the projection for the sampling coordinates:

coordinates(lonlat) <- ~lon + lat
proj4string(lonlat) <- proj4string(income_merged)

# Extract the census data to the sampling coordinates:

llData <- over(lonlat, income_merged) %>%
  cbind(lonlat@data) %>%
  cbind(lonlat@coords) %>%
  tbl_df %>%
  select(site, lon, lat, GEOID, can, imp, medianIncome) %>%
  dplyr::rename(fips = GEOID)

# Take a look:

hist(llData$medianIncome)

# Write to file:

write.csv(llData, 'llData.csv', row.names = FALSE)

