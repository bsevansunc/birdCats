library(RCurl)
library(stringr)
library(tidyverse)

x <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/captureTable.csv'
  )

y <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/nnPowerAnalysis/master/nnData/measurementTable.csv'
  )

z <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/birdCats/master/birdCats/data/catDataTransect.csv'
)

sites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1','GERYERIMD1', 'MISSEDDC1')


addBand <- read.csv(text = x) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(
    str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
    site %in% sites,
    str_detect(birdID, '[0-9]{3}-[0-9]{5}')|
      str_detect(birdID, '[0-9]{4}-[0-9]{5}'),
    typeCapture== 'Band') %>%
  mutate(birdID = as.character(birdID)) %>%
  arrange(birdID) %>%
  select(birdID, site, date, typeCapture, age, sex) %>%
  distinct

addRecap <- read.csv(text = x) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  mutate(birdID = as.character(birdID)) %>%
  filter(
    str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
    site %in% sites,
    str_detect(birdID, '[0-9]{3}-[0-9]{5}')|
      str_detect(birdID, '[0-9]{4}-[0-9]{5}'),
    typeCapture == 'Recap') %>%
  arrange(birdID) %>%
  select(birdID, site, date, typeCapture) %>%
  distinct

birdIDs <- distinct(select(addBand,birdID))[,1]

addMeas <- read.csv(text = y) %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  mutate(birdID = as.character(birdID)) %>%
  filter(str_detect(date, '[0-9]{4}-[0-9]{2}-[0-9]{2}'),
         site %in% sites,
         str_detect(birdID, '[0-9]{3}-[0-9]{5}')|
           str_detect(birdID, '[0-9]{4}-[0-9]{5}'),
         !is.na(mass) & !is.na(wing),
         birdID %in% birdIDs) %>%
  arrange(birdID,date) %>%
  select(birdID, mass, wing) %>%
  distinct(birdID,.keep_all = TRUE)

addBandMeas <- addBand %>%
  left_join(addMeas, by='birdID') %>%
  distinct

write.csv(addRecap, 'data/addRecap.csv')
write.csv(addBandMeas, 'data/addBand.csv')



