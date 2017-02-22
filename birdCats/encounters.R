library('tidyverse')
library('stringr')
library('lubridate')

# --------------------------------------------------------*
# ----------- Read in and subset the data ----------------
# --------------------------------------------------------*

# Desired species

species <- c('AMRO', 'CACH', 'CARW', 'GRCA',
             'HOWR', 'NOCA', 'NOMO', 'SOSP')

# Sites to remove from analyses:

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1','GERYERIMD1')


# Read in banding and recap data, filter to observations with band number, 
# band encounters, and birds banded before 2016. Subset columns to
# bandNumber, site, year, species, age, and sex:

band <- read_csv('encBand.csv') %>%
  arrange(bandNumber, date) %>%
  mutate(age = str_replace_all(age, 'UNK', 'U'),
         sex = str_replace_all(sex, 'UNK', 'U'),
         year = year(date)) %>%
  filter(
    str_detect(bandNumber, '[0-9]{3}-[0-9]{5}')|
      str_detect(bandNumber, '[0-9]{4}-[0-9]{5}'),
    speciesEnc %in% species,
    encounterType == 'Band',
    year != 2016,
    !site %in% removeSites
    ) %>%
  select(bandNumber, site, year, speciesEnc, age, sex) %>%
  distinct

# Read in banding and recap data, filter observations without a band number, 
# recap encounters only, and subset columns to bandNumber and year:

recap <- read_csv('encBand.csv') %>%
  arrange(bandNumber, date) %>%
  mutate(year = year(date)) %>%
  filter(bandNumber %in% band$bandNumber) %>%
  select(bandNumber, year) %>%
  distinct

# Read in participant resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

partResight <- read_csv('encPart.csv') %>%
  filter(birdID %in% band$bandNumber) %>%
  # The below ensures that resights are roughly 1 year after capture:
  mutate(yearResight = ifelse(month(dateResight) < 4,
                              yearResight -1, yearResight)) %>%
  arrange(birdID, yearResight) %>%
  transmute(bandNumber = birdID, year = yearResight) %>%
  distinct

# Read in technician resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

techResight <- read_csv('encTech.csv') %>%
  separate(visitID, c('site', 'date'), '\\_') %>%
  filter(birdID %in% bandRecap$bandNumber) %>%
  transmute(bandNumber = birdID,
            year = year(date)) %>%
  arrange(bandNumber, year) %>%
  distinct

# Read in K. Bennet resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

kbResight <- read_csv('encKB.csv') %>%
  filter(band %in% bandRecap$bandNumber) %>%
  transmute(bandNumber = band,
            year = year(date)) %>%
  arrange(bandNumber, year) %>%
  distinct
  
# Combine all encounter types to make capture history frame:

captureHistories <- bind_rows(
  band %>% select(bandNumber, year),
  recap, partResight, techResight, kbResight) %>%
  arrange(bandNumber, year) %>%
  mutate(enc = 1) %>%
  distinct %>%
  spread(key = year, value = enc, fill = 0) %>%
  unite(ch, -bandNumber, sep = '')

# Make data frame that includes capture histories, groups, and covariates:

birdData <- left_join(bandRecap,captureHistories, by = 'bandNumber') %>%
  rename(bandYear = year)







