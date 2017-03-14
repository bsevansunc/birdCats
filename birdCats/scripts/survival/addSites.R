library(RCurl)
library(stringr)
library(tidyverse)

u <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/birdTable.csv'
)

x <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/captureTable.csv'
  )

y <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/measurementTable.csv'
  )

z <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/birdCats/master/birdCats/data/catDataTransect.csv'
)

v <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/resightPartTable.csv'
)

w <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/resightTechTable.csv'
)

sites <- read.csv('data/catDataTransect.csv')$site


getBand <- read.csv(text = x) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(
    str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
    site %in% sites,
    str_detect(bandNumber, '[0-9]{3}-[0-9]{5}')|
      str_detect(bandNumber, '[0-9]{4}-[0-9]{5}'),
    typeCapture== 'Band') %>%
  mutate(year = year(date)) %>%
  filter(year != 2016) %>%
  arrange(bandNumber) %>%
  select(bandNumber, site, date, age, sex) %>%
  distinct

getRecap <- read.csv(text = x) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(
    str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
    site %in% sites,
    str_detect(bandNumber, '[0-9]{3}-[0-9]{5}')|
      str_detect(bandNumber, '[0-9]{4}-[0-9]{5}'),
    typeCapture == 'Recap') %>%
  arrange(bandNumber) %>%
  select(bandNumber, site, date) %>%
  distinct

getMeas <- read.csv(text = y) %>%
  mutate(bandNumber = birdID) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
         site %in% sites,
         str_detect(bandNumber, '[0-9]{3}-[0-9]{5}')|
           str_detect(bandNumber, '[0-9]{4}-[0-9]{5}'),
         !is.na(mass) & !is.na(wing),
         bandNumber %in% getBand$bandNumber) %>%
  arrange(bandNumber,date) %>%
  select(bandNumber, mass, wing) %>%
  distinct(bandNumber,.keep_all = TRUE)


spp <- read.csv(text=u) %>%
  select(bandNumber,species) %>%
  filter(bandNumber %in% getBand$bandNumber) %>%
  distinct


bandMeas <- getBand %>%
  left_join(getMeas, by='bandNumber') %>%
  left_join(spp, by='bandNumber') %>%
  distinct


getPart <- read.csv(text = v) %>%
  filter(siteID %in% sites) %>%
  select(siteID, birdID, dateResight)

getTech <- read.csv(text = w) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(site %in% sites) %>%
  select(site,date,birdID)


write.csv(getRecap, 'data/encRecap.csv',row.names = FALSE)
write.csv(bandMeas, 'data/encBand.csv',row.names = FALSE)
write.csv(getPart, 'data/encPart.csv',row.names = FALSE)
write.csv(getTech, 'data/encTech.csv',row.names = FALSE)



