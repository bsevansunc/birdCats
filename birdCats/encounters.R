library('dplyr')
library('tidyr')

setwd('C:/Users/kbenn/Documents/GitHub/birdCats/birdCats')



# --------------------------------------------------------*
# ----------- Read in and subset the data ----------------
# --------------------------------------------------------*

# Desired species

species <- c('AMRO', 'CACH', 'CARW', 'GRCA', 'HOWR', 'NOCA', 'NOMO', 'SOSP')



# Read in the data, subset banding data to site, date, species, band, age, and sex,
# and subset the resight data to band number and date

band <- read.csv('encBand.csv') %>%
  select(bandNumber, site, date, speciesEnc, age, sex) %>%
  arrange(bandNumber, date) %>%
  lapply(gsub, pattern = 'UNK', replacement = 'U') %>%
  data.frame() %>%
  mutate(bandNumber = as.character(bandNumber)) %>%
  filter(bandNumber != 't-estBn' & bandNumber != 'tes-tBand' & !is.na(bandNumber)) %>%
  filter(speciesEnc %in% species) %>%
  unique()

part <- read.csv('encPart.csv') %>%
  select(birdID, yearResight) %>%
  arrange(birdID, yearResight) %>%
  rename(bandNumber = birdID, year = yearResight) %>%
  mutate(bandNumber = as.character(bandNumber)) %>%
  filter(!is.na(bandNumber)) %>%
  unique()

tech <- read.csv('encTech.csv') %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  select(birdID, date) %>%
  arrange(birdID, date) %>%
  rename(bandNumber = birdID) %>%
  mutate(bandNumber = as.character(bandNumber)) %>%
  filter(bandNumber != '-') %>%
  unique()

kb <- read.csv('encKB.csv') %>%
  filter(site != 'OLONMARDC1' & site != 'WOLFKARDC1' & 
         site != 'WOLFAMYDC1' & site != 'GERYERIMD1') %>%
  select(band, date) %>%
  arrange(band, date) %>%
  rename(bandNumber = band) %>%
  mutate(bandNumber = as.character(bandNumber)) %>%
  filter(bandNumber != '?') %>%
  unique()


# Create column 'year' from band date

band$date <- as.Date(band$date)
band$year <- as.integer(format(band$date, '%Y'))
band$date <- NULL

tech$date <- as.Date(tech$date)
tech$year <- as.integer(format(tech$date, '%Y'))
tech$date <- NULL

kb$date <- as.Date(kb$date)
kb$year <- as.integer(format(kb$date, '%Y'))
kb$date <- NULL


# Combine

all <- bind_rows(list(band, tech, part, kb)) %>%
  arrange(bandNumber, year) %>%
  select(bandNumber, year) %>%
  mutate(enc = 1) %>%
  unique()



# --------------------------------------------------------*
# -------------- Create encounter history ----------------
# --------------------------------------------------------*


# Spread to encounter history by individual, 
# then get frequencies of each encounter history

enc <- spread(all, year, enc, fill = 0) %>%
  unite(ch, -bandNumber, sep = '') %>%
  arrange(ch) %>%
  select(ch) %>%
  table() %>%
  data.frame() %>%
  setNames(c('ch', 'freq'))


# Write to .csv file

write.csv(enc, 'encounterHistory.csv', row.names = FALSE)




