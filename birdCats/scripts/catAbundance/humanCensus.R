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
                  table.number = "B19013", col.names = "pretty")

# Make into a data frame:

income_df <- data.frame(
  GEOID = paste0(str_pad(income@geography$state, 2, "left", pad="0"),
                 str_pad(income@geography$county, 3, "left", pad="0"),
                 str_pad(income@geography$tract, 6, "left", pad="0")),
  medianIncome = as.numeric(income@estimate)
  )

population <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                        table.number = "B01003", col.names = "pretty")
  
population_df <- data.frame(
  GEOID = paste0(str_pad(population@geography$state, 2, "left", pad="0"),
                 str_pad(population@geography$county, 3, "left", pad="0"),
                 str_pad(population@geography$tract, 6, "left", pad="0")),
  population = as.numeric(population@estimate)
) %>%
  left_join(
    data.frame(
      GEOID = tractsAll$GEOID,
      area = tractsAll$ALAND
    )
  ) %>%
  mutate(hDensity = population/area)

age <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                 table.number = "B01002", col.names = "pretty")

age_df <- data.frame(
  GEOID = paste0(str_pad(age@geography$state, 2, "left", pad="0"),
                 str_pad(age@geography$county, 3, "left", pad="0"),
                 str_pad(age@geography$tract, 6, "left", pad="0")),
  age = as.numeric(age@estimate)
)

edu <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                 table.number = "B15003", col.names = "pretty")

edu_df <- data.frame(
  GEOID = paste0(str_pad(edu@geography$state, 2, "left", pad="0"),
                 str_pad(edu@geography$county, 3, "left", pad="0"),
                 str_pad(edu@geography$tract, 6, "left", pad="0")),
  eduHS = (edu@estimate[,17])/edu@estimate[,1],
  eduC = (edu@estimate[,22])/edu@estimate[,1]
)


marred <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                 table.number = "B12001", col.names = "pretty")

marred_df <- data.frame(
  GEOID = paste0(str_pad(marred@geography$state, 2, "left", pad="0"),
                 str_pad(marred@geography$county, 3, "left", pad="0"),
                 str_pad(marred@geography$tract, 6, "left", pad="0")),
  marred = (marred@estimate[,5] + marred@estimate[,13])/marred@estimate[,1]
)


# Merge census data with geographic dataset:

income_merged<- geo_join(tractsAll, income_df, "GEOID", "GEOID")
population_merged<- geo_join(tractsAll, population_df, "GEOID", "GEOID")
age_merged<- geo_join(tractsAll, age_df, "GEOID", "GEOID")
edu_merged<- geo_join(tractsAll, edu_df, "GEOID", "GEOID")
marred_merged<- geo_join(tractsAll, marred_df, "GEOID", "GEOID")





# ---------------------------------------------------------------------------------*
# ---- Plot census data ----
# ---------------------------------------------------------------------------------*

# Popup:

popup <- paste0("GEOID: ", income_merged$GEOID, "<br>", "Median household income: ")

# Color palette:

pal <- colorNumeric(
  palette = "YlGnBu",
  domain = population_merged$hDensity
  # domain = income_merged$medianIncome
)

map <-leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = population_merged, #income_merged, 
              fillColor = ~pal(hDensity), #~pal(medianIncome), 
              color = "#b2aeae", # you need to use hex colors
              fillOpacity = 0.7, 
              weight = 1, 
              smoothFactor = 0.2,
              popup = popup) %>%
  addLegend(pal = pal, 
            values = population_merged$hDensity, #income_merged$medianIncome, 
            position = "bottomright", 
            title = "Median household<br>income ($)")#,
            # labFormat = labelFormat(prefix = "")) 
map

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
  filter(!is.na(cat)) %>%
  select(site) %>%
  distinct %>%
  arrange(site) %>%
  .$site

# Transform sampling coordinates from UTM to lonlat and filter to cat sampling points:

lonlat <- spTransform(utm, CRS('+proj=longlat')) %>%
  data.frame %>%
  dplyr::rename(lon = coords.x1, lat = coords.x2) %>%
  cbind(lc100) %>%
  select(site, can, imp, lon, lat) %>%
  filter(site %in% catSites) %>%
  bind_rows(
    data.frame(
      site = c('GERYERIMD1', 'MISSEDDC1', 'OBRICHRMD1',
               'OLONMARDC1', 'WOLFAMYDC1', 'WOLFKARDC1'),
      lon = c(-77.118616, -76.996958, -76.924433,
              -77.083900, -77.043167,  -77.061800),
      lat = c(38.999316, 38.883588,39.013183,
              38.949000,38.918567, 38.941583)
    )
  ) %>%
  data.frame

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

siteCensusInfo <- function(mergeFile, variable){
  over(lonlat, mergeFile) %>%
    cbind(lonlat@data) %>%
    cbind(lonlat@coords) %>%
    tbl_df %>%
    select_('site', variable)
}

mergeList <- list(population_merged, age_merged, edu_merged, edu_merged, marred_merged)
variables <- c('hDensity', 'age', 'eduHS','eduC', 'marred')

outList  <- vector('list', length = length(mergeList))

for(i in 1:length(mergeList)){
  outList[[i]] <- siteCensusInfo(mergeList[[i]], variables[i])
}

llDataCensus <- llData %>%
  left_join(outList[[1]]) %>%
  left_join(outList[[2]]) %>%
  left_join(outList[[3]]) %>%
  left_join(outList[[4]]) %>%
  left_join(outList[[5]]) %>%
  mutate(can = ifelse(site == 'GERYERIMD1', 21.482, can),
         imp = ifelse(site == 'GERYERIMD1', 32.225, imp),
         can = ifelse(site == 'MISSEDDC1', 2.36, can),
         imp = ifelse(site == 'MISSEDDC1', 61.429, imp),
         can = ifelse(site == 'OBRICHRMD1', 29.381, can),
         imp = ifelse(site == 'OBRICHRMD1', 26.456, imp),
         can = ifelse(site == 'OLONMARDC1', 26.327, can),
         imp = ifelse(site == 'OLONMARDC1', 29.959, imp),
         can = ifelse(site == 'WOLFAMYDC1', 0, can),
         imp = ifelse(site == 'WOLFAMYDC1', 81.033, imp),
         can = ifelse(site == 'WOLFKARDC1', 5.22, can),
         imp = ifelse(site == 'WOLFKARDC1', 48.631, imp)
         )

  
# ptsCanImp <- read.csv('pts.can.imp.6.6.csv', stringsAsFactors = FALSE) %>%
#   select(SITE, can100, imp100) %>%
#   rename(site = SITE)
# 
# llDataCensus %>%
#   left_join(ptsCanImp, by = 'site') %>% View


# Take a look:

hist(llData$medianIncome)

# Write to file:

write.csv(llDataCensus, 'covariateData.csv', row.names = FALSE)

