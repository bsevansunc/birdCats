library(dplyr) ; library(ggplot2) ; library(ggmap)

kmlPath <- "bandingStation_upperZoo.kml"

kmlToDf <- function(kmlPath){
  require(plyr) ; require(maptools)
  getKMLcoordinates(kmlfile=kmlPath, ignoreAltitude=T) %>%
    ldply %>%
    rename(lon = V1, lat = V2)
}

pointFrame <- bind_rows(
  kmlToDf('bandingStation_upperZoo.kml') %>%
    mutate(class = 'bandingStation',
           location = 'upper'),
  kmlToDf('bandingStation_lowerZoo.kml') %>%
    mutate(class = 'bandingStation',
           location = 'lower'),
  kmlToDf('netLanes_upperZoo.kml') %>%
    mutate(class = 'netLanes',
           location = 'upper'),
  kmlToDf('netLanes_lowerZoo.kml') %>%
    mutate(class = 'netLanes',
           location = 'lower')
  )

ggplot(pointFrame, aes(x = lat, y = lon, color = location, shape = class)) + geom_point()

gMap <- get_googlemap(center = c(-77.049, 38.93), zoom = 16, scale = 2)

ggmap(gMap)
plot(gMap)


