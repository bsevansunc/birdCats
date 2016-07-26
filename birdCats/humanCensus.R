library(tigris)
library(acs)
library(stringr) # to pad fips codes
library(leaflet)

api.key.install(key = 'ec7729865a6a6bdf5d1ce20999a4788da99b6c5e')

countiesMD <- c(3,5,510, 31, 33)
countiesVA <- c(600, 13, 510, 113, 179)
countiesDC <- 11001


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

latlong2fips <- function(latitude, longitude) {
  url <- "http://data.fcc.gov/api/block/find?format=json&latitude=%f&longitude=%f"
  url <- sprintf(url, latitude, longitude)
  json <- RCurl::getURL(url)
  json <- RJSONIO::fromJSON(json)
  as.character(json$County['FIPS'])
}




