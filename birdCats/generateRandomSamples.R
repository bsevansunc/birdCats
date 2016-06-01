# Set-up:

# Working directory, Brian office: C:/Users/Brian/Desktop/gits/birdCats
# Working directory, Kevin's laptop: C:/Users/Kevin/Documents/GitHub/birdCats/birdCats

library(RCurl)
library(dplyr)
library(tidyr)

# Read sites for study:

gitURL <- getURL(
  'https://raw.githubusercontent.com/bsevansunc/birdCats/master/birdCats/catSites.csv'
)

# Get sites:

catSites <- read.csv(text = gitURL, stringsAsFactors = FALSE) %>%
  tbl_df %>%
  .$site

# Cat frame with month (1:3) and sample order (1:53, per month):


catFrame <- data.frame(month = rep(1:3, each = length(unique(catSites))),
           sampleOrder = rep(1:53,  3),
           site = NA) %>%
  tbl_df

# Sites for month 1 are a random sample of catSites:

catFrame[1:53,'site'] <- sample(catSites, 53)

# Get vector of the last 14 cat sites in month 1:

endCats1 <- catFrame %>%
  filter(month == 1, sampleOrder >= 40) %>%
  .$site

# Generate random sample of sites for first 14 sites in month 2 that does not include endCats1

catFrame[54:67,'site'] <- sample(
  catSites[!catSites %in% endCats1], 14
)

# For remaining sites in month 2, select from vector of sites that do not include the first 14 records of this month:

catFrame[68:106,'site'] <- sample(
  catSites[!catSites %in% catFrame[54:67,]$site], 39 
)

# Get vector of the last 14 cat sites in month 2:

endCats2 <- catFrame %>%
  filter(month == 2, sampleOrder >= 40) %>%
  .$site

# Generate random sample of sites for first 14 sites in month 3 that does not include endCats2

catFrame[107:120,'site'] <- sample(
  catSites[!catSites %in% endCats2], 14
)

# For remaining sites in month 3, select from vector of sites that do not include the first 14 records of this month:

catFrame[121:nrow(catFrame),'site'] <- sample(
  catSites[!catSites %in% catFrame[107:120,]$site], 39 
)

# Write file:

write.csv(catFrame, 'samplingOrder2.csv', row.names = FALSE)


