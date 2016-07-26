library(tigris)
library(acs)
library(stringr) # to pad fips codes
library(leaflet)
library(rgdal)

api.key.install(key = 'ec7729865a6a6bdf5d1ce20999a4788da99b6c5e')




tractsMD <- tracts(state = 'MD')
tractsVA <- tracts(state = 'VA')
tractsDC <- tracts(state = 'DC')

tractsAll <- rbind_tigris(tractsMD, tractsVA, tractsDC)

leaflet(tractsAll) %>%
  addTiles() %>%
  addPolygons(popup = ~NAME)

geo <- geo.make(state = c('DC', 'VA', 'MD'), county = '*', tract = '*')

income <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                  table.number = "B19001", col.names = "pretty")

income1 <- acs.fetch(endyear = 2012, span = 5, geography = geo,
                  table.number = "B19013", col.names = "pretty")

names(attributes(income))

income_df <- data.frame(paste0(str_pad(income@geography$state, 2, "left", pad="0"), 
                               str_pad(income@geography$county, 3, "left", pad="0"), 
                               str_pad(income@geography$tract, 6, "left", pad="0")),
                        as.numeric(income1@estimate),
#                         income@estimate[,c("Household Income: Total:",
#                                            "Household Income: $200,000 or more")], 
                        stringsAsFactors = FALSE)
# 
# income_df <- select(income_df, 1:3)
# rownames(income_df)<-1:nrow(income_df)
# names(income_df)<-c("GEOID", "total", "over_200")
# income_df$percent <- 100*(income_df$over_200/income_df$total)

names(income_df) <- c('GEOID', 'medianIncome')

income_merged<- geo_join(tractsAll, income_df, "GEOID", "GEOID")
# there are some tracts with no land that we should exclude
income_merged <- income_merged[income_merged$ALAND>0,]

popup <- paste0("GEOID: ", income_merged$GEOID, "<br>", "Median household income: ")
# popup <- paste0("GEOID: ", income_merged$GEOID, "<br>", "Percent of Households above $200k: ", round(income_merged$percent,2))
pal <- colorNumeric(
  palette = "YlGnBu",
  domain = income_merged$medianIncome
  # domain = income_merged$percent
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

# Add to census data to cat sampling points

lc100 <- read.csv('lc100.csv')

utm <- SpatialPoints(cbind(lc100$X,lc100$Y), proj4string=CRS("+proj=utm +zone=18"))

catSites <- read.csv('catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT)) %>%
  select(SITE) %>%
  distinct %>%
  arrange(SITE) %>%
  .$SITE

lonlat <- spTransform(utm, CRS('+proj=longlat')) %>%
  data.frame %>%
  dplyr::rename(lon = coords.x1, lat = coords.x2) %>%
  cbind(lc100) %>%
  select(site, can, imp, lon, lat) %>%
  filter(site %in% catSites)

coordinates(lonlat) <- ~lon + lat
proj4string(lonlat) <- proj4string(income_merged)

llData <- over(lonlat, income_merged) %>%
  cbind(lonlat@data) %>%
  tbl_df %>%
  select(site, GEOID, can, imp, medianIncome) %>%
  dplyr::rename(fips = GEOID)

hist(llData$medianIncome)
  




