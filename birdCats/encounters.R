library('dplyr')
library('tidyr')

setwd('C:/Users/kbenn/Documents/GitHub/birdCats/birdCats')



# --------------------------------------------------------*
# ---- Read in the data and preliminary manipulations ----
# --------------------------------------------------------*

# Read in the data

band <- read.csv('encBand.csv')
part <- read.csv('encPart.csv')


# Create column 'year' from band date

band$date <- as.Date(band$date)
band$year <- as.integer(format(band$date, '%Y'))


# Arrange both dataframes by band number and year, and
# make band numbers character strings

band <- arrange(band, bandNumber, year)
part <- arrange(part, birdID, yearResight)

band$bandNumber <- as.character(band$bandNumber)
part$birdID <- as.character(part$birdID)



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



# --------------------------------------------------------*
# --------- Subset to desired data and combine -----------
# --------------------------------------------------------*

# Desired species

species <- c('AMRO', 'CACH', 'CARW', 'GRCA', 'HOWR', 'NOCA', 'NOMO', 'SOSP')


# Subset to desired columns

a <- band %>%
  filter(speciesEnc %in% species) %>%
  select(bandNumber, site, year)

b <- part %>%
  select(birdID, siteID, yearResight) %>%
  rename(bandNumber = birdID, site = siteID, year = yearResight)


# Keep sites as a factor

levels(b$site) <- levels(a$site)


# Combine

c <- bind_rows(a, b) %>%
  arrange(bandNumber, year) %>%
  unique() %>%
  mutate(enc = 1)



# --------------------------------------------------------*
# -------------- Create encounter history ----------------
# --------------------------------------------------------*

enc <- spread(c, year, enc, fill = 0)

write.csv(enc, 'encounterHistory.csv', row.names = FALSE)

