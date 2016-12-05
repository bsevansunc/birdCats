library(shiny)
library(dplyr)
library(tidyr)
library(leaflet)

function(input, output) {
  # Set color pallete:
  
  pal <- colorFactor(c("navy", "red"), domain = c("net", "banding"))
  
  # Generate map:
  map <- leaflet(data = netPoints) %>%
    addTiles() %>%
    addCircleMarkers(~long, ~lat,
                     radius = ~ifelse(ptType == 'net', 6, 10),
                     col = ~pal(ptType),
                     stroke = FALSE,
                     fillOpacity = 0.8,
                     popup = ~as.character(ptName))
  output$zooMap <- renderLeaflet(map)
}