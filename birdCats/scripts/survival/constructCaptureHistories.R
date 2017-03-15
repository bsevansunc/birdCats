library('tidyverse')
library('stringr')
library('lubridate')
library('scales')

# --------------------------------------------------------*
# ----------- Read in and subset the data ----------------
# --------------------------------------------------------*

# Desired species

studySpecies <- c('AMRO', 'CACH', 'CARW', 'GRCA',
             'HOWR', 'NOCA', 'NOMO', 'SOSP')


# Sites to remove from analyses:

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1','GERYERIMD1', 'MISSEDDC1')


# Read in banding data

band <- read_csv('data/encBand.csv') %>%
  arrange(bandNumber, date) %>%
  mutate(age = str_replace_all(age, 'UNK', 'U'),
         sex = str_replace_all(sex, 'UNK', 'U'),
         year = year(date)) %>%
  filter(species %in% studySpecies) %>%
  distinct


# Read in recap data. Subset columns to bandNumber, site, and year:

recap <- read_csv('data/encRecap.csv') %>%
  arrange(bandNumber, date) %>%
  mutate(year = year(date)) %>%
  select(bandNumber, year, site) %>%
  distinct


# Read in participant resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

partResight <- read_csv('data/encPart.csv') %>%
  filter(birdID %in% band$bandNumber) %>%
  # The below ensures that resights are roughly 1 year after capture:
  mutate(
    yearResight = year(dateResight),
    yearResight = ifelse(month(dateResight) < 4,
                         yearResight -1, yearResight)) %>%
  arrange(birdID, yearResight) %>%
  transmute(site = siteID, bandNumber = birdID, year = yearResight) %>%
  distinct


# Read in technician resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

techResight <- read_csv('data/encTech.csv') %>%
  filter(birdID %in% band$bandNumber) %>%
  transmute(site = site, bandNumber = birdID, year = year(date)) %>%
  arrange(bandNumber, year) %>%
  distinct


# Read in K. Bennett resights, filter to band numbers in band data, and
# subset columns to bandNumber and year of resight:

kbResight <- read_csv('data/encKB.csv') %>%
  filter(band %in% band$bandNumber) %>%
  transmute(bandNumber = band, year = year(date), site = site) %>%
  arrange(bandNumber, year) %>%
  distinct



# ---------------------------------------------------------------------*
# ---- Create body condition variable based on Peig and Green 2009 ----
# ---------------------------------------------------------------------*

# Function that creates a dataframe of ID, mass, wing, sex per species
sppDf <- function(df, sp){
  df %>% filter(species == sp) %>% select(bandNumber, mass, wing, sex)
}


# Make the dataframes

amro <- sppDf(band, 'AMRO')
cach <- sppDf(band, 'CACH')
carw <- sppDf(band, 'CARW')
grca <- sppDf(band, 'GRCA')
howr <- sppDf(band, 'HOWR')
noca <- sppDf(band, 'NOCA')
nomo <- sppDf(band, 'NOMO')
sosp <- sppDf(band, 'SOSP')


# Remove outliers in mass and wing data

rmOutliers <- function(df){
  df$mw <- round(df$mass/df$wing, digits=8)
  repeat{
    if (grubbs.test(df$mw)$p.value > 0.05) break
    out <- round(
      as.numeric(gsub("[^0-9.]","",grubbs.test(df$mw)$alternative)),digits=8)
    df <- df %>% filter(mw != out | is.na(mw))
  }
  return(df %>% select(-mw))
}

amro <- rmOutliers(amro)
cach <- rmOutliers(cach)
carw <- rmOutliers(carw)
grca <- rmOutliers(grca)
howr <- rmOutliers(howr)
noca <- rmOutliers(noca)
nomo <- rmOutliers(nomo)
sosp <- rmOutliers(sosp)



#Function that calculates the body condition index

bciFun <- function(species){
  df1 <- species
  L0 <- mean(df1$wing, na.rm = TRUE)
  Mi <- df1$mass
  Li <- df1$wing
  bsma <- summary(lm(log(Mi)~log(Li)))$coefficients[2,1]/cor(
    log(Mi),log(Li), use = 'complete.obs')

  Mi*((L0/Li)^bsma)
}


# Calculate BCI for each species, scaled from 0 to 1

amro$bci <- rescale(bciFun(amro))
cach$bci <- rescale(bciFun(cach))
carw$bci <- rescale(bciFun(carw))
grca$bci <- rescale(bciFun(grca))
howr$bci <- rescale(bciFun(howr))
noca$bci <- rescale(bciFun(noca))
nomo$bci <- rescale(bciFun(nomo))
sosp$bci <- rescale(bciFun(sosp))


#Combine into one dataframe

bciDf <- bind_rows(
  amro %>% select(bandNumber, bci),
  cach %>% select(bandNumber, bci),
  carw %>% select(bandNumber, bci),
  grca %>% select(bandNumber, bci),
  howr %>% select(bandNumber, bci),
  noca %>% select(bandNumber, bci),
  nomo %>% select(bandNumber, bci),
  sosp %>% select(bandNumber, bci))


# Add to band dataframe with other covariates

band <- left_join(band, bciDf, by = 'bandNumber')


# -------------------------------------------------------*
# ----------- Add urbanization variables -----------
# -------------------------------------------------------*

# Read in data

urban <- read_csv('data/covariateData.csv') %>%
  select(site, imp, can) %>%
  filter(site %in% band$site)


# Add to band dataframe with other covariates

band <- left_join(band, urban, by = 'site')


# ------------------------------------------------------------------*
# ---------------- Add cat variables -------------------------------
# ------------------------------------------------------------------*

# Read in data

catAbund <- read.csv('data/catAbund.csv')


# Add to band dataframe with other covariates

band <- left_join(band, catAbund, by = 'site')



# --------------------------------------------------------------------*
# ------------ Identify sites as active or inactive ------------------
# --------------------------------------------------------------------*

# Create a dataframe of sites and years in which 
# banding/resighting occurred at each

siteHistory <- bind_rows(band %>% select(site, year),
                         recap %>% select(site, year),
                         partResight %>% select(site, year),
                         techResight %>% select(site, year),
                         kbResight %>% select(site, year)) %>%
  arrange(site, year) %>%
  distinct


# Create a dataframe of bird ID and corresponding site

idSite <- band %>% select(bandNumber, site)


# Create a wide dataframe with '1' indicating an active site,
# '0' indicating inactive

yearActive <- left_join(idSite, siteHistory, by = 'site') %>%
  mutate(active = '1') %>%
  spread(key = year, value = active, fill = 0) %>%
  select(-site)

colnames(yearActive) <- c('bandNumber', paste('active', 1:17, sep = ''))


# ---------------------------------------------------------*
# ----------- Construct capture histories -----------------
# ---------------------------------------------------------*


# Combine all encounter types to make capture history frame:

captureHistories <- bind_rows(
  band %>% select(bandNumber, year),
  recap %>% select(-site),
  partResight %>% select(-site),
  techResight %>% select(-site), 
  kbResight %>% select(-site)) %>%
  arrange(bandNumber, year) %>%
  mutate(enc = 1) %>%
  distinct %>%
  spread(key = year, value = enc, fill = 0) %>%
  unite(ch, -bandNumber, sep = '')



# Make data frame that includes capture histories, groups, and covariates

birdData <- left_join(band,captureHistories, by = 'bandNumber') %>%
  left_join(yearActive, by = 'bandNumber') %>%
  rename(bandDate = date, bandYear = year)

# Writes to file for survival analysis

write.csv(birdData, 'data/ch.csv', row.names = FALSE)
