library('dplyr')
library('tidyr')

setwd('C:/Users/kbenn/Documents/GitHub/birdCats/birdCats')



# --------------------------------------------------------*
# ---- Read in the data and preliminary manipulations ----
# --------------------------------------------------------*

# Read in the data

band <- read.csv('encBand.csv')
part <- read.csv('encPart.csv')
tech <- read.csv('encTech.csv')
kb <- read.csv('encKB.csv')


# New variables site and date from visitID in tech

tech <- separate(tech, visitID, c('site', 'date'), '\\_')
tech$site <- as.factor(tech$site)


# Create column 'year' from band date

band$date <- as.Date(band$date)
band$year <- as.integer(format(band$date, '%Y'))

tech$date <- as.Date(tech$date)
tech$year <- as.integer(format(tech$date, '%Y'))

kb$date <- as.Date(kb$date)
kb$year <- as.integer(format(kb$date, '%Y'))


# Arrange dataframes by band number and year, and
# make band numbers character strings

band <- arrange(band, bandNumber, year)
part <- arrange(part, birdID, yearResight)
tech <- arrange(tech, birdID, year)
kb <- arrange(kb, band, year)

band$bandNumber <- as.character(band$bandNumber)
part$birdID <- as.character(part$birdID)
tech$birdID <- as.character(tech$birdID)
kb$band <- as.character(kb$band)



# --------------------------------------------------------*
# -------- Fix errors and clean up the dataset -----------
# --------------------------------------------------------*

# Band entered incorrectly

band[2745,6] <- '2711-77659'
part[617,4] <- '2641-63884'


# Remove NAs and tests

band <- filter(
  band, bandNumber != '-' & bandNumber != 't-estBn' & bandNumber != 'tes-tBand')

part <- filter(part, !is.na(birdID))

tech <- filter(tech, birdID != '-')

kb <- filter(kb, band != '?')



# --------------------------------------------------------*
# --------- Subset to desired data and combine -----------
# --------------------------------------------------------*

# Desired species

species <- c('AMRO', 'CACH', 'CARW', 'GRCA', 'HOWR', 'NOCA', 'NOMO', 'SOSP')


# Subset to desired columns

a <- band %>%
  filter(speciesEnc %in% species) %>%
  select(bandNumber, year)

b <- part %>%
  select(birdID, yearResight) %>%
  rename(bandNumber = birdID, year = yearResight)

c <- tech %>%
  select(birdID, year) %>%
  rename(bandNumber = birdID)

d <- kb %>%
  select(band, year) %>%
  rename(bandNumber = band)


# Combine

e <- bind_rows(list(a, b, c, d)) %>%
  arrange(bandNumber, year) %>%
  unique() %>%
  mutate(enc = 1)



# --------------------------------------------------------*
# -------------- Create encounter history ----------------
# --------------------------------------------------------*

# Encounter history by individual

encInd <- spread(e, year, enc, fill = 0) %>%
  unite(ch, -bandNumber, sep = '') %>%
  arrange(ch) %>%
  select(ch)


# Frequencies of encounter histories

enc <- as.data.frame(table(encInd)) %>%
  rename(ch = encInd, freq = Freq)


# Write to .csv file

write.csv(enc, 'encounterHistory.csv', row.names = FALSE)




