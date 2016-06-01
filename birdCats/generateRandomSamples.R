# Set-up:

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

# Generate site data frame, with no week 4 to week 1 transitions:

t1 <- data.frame(site = rep(catSites,3)) %>%
  arrange(site) %>%
  # For each site, assign sample months 1,2,3:
  mutate(month = rep(1:3, length(unique(site)))) %>%
  # Grouping operation allows us to sample within site and month:
  group_by(site, month) %>%
  # Add a random "week" sample of 1-4:
  mutate(week = sample(1:4, 1)) %>%
  # Bookkeeping for better visualization
  ungroup %>%
  arrange(site) 

# Empty list, with 1 list item per site in which to store data:
tList <- vector('list', length(catSites))

# For each site ...
for(i in 1:length(tList)){
  # Data frame of one site, arranged by month:
  siteItem <- filter(t1, site == catSites[[i]]) %>%
    arrange(month)
  # For each month sample except the first  ...
  for(j in 2:3){
    # If the sample week is 1 and the previous week is 4, resample from 2:4
    if(siteItem[j,3] == 1 & siteItem[j-1,3] == 4){
      siteItem[j, 3] <- sample(2:4, 1)
    }
  }
  # Output list:
  tList[[i]] <- siteItem
}

# Bind the fixed list back into a data frame:

t1 <- bind_rows(tList)

# Make an empty list of sample events by month:

t1List <- vector('list', length = 3)

# For each month ...

for(i in 1:3){
  # Filter frame to that month
  t1Month <- t1 %>%
    filter(month == i)
  # Get a vector of "weeks" within that month:
  weekList <- vector('list', length = length(unique(t1Month$week)))
  # For each week in month ...
  for(j in 1:length(weekList)){
    # Filter frame to that week:
    t1Week <- filter(t1Month, week == j)
    # Add ordered sample events within that week:
    t1Week$event <- sample(1:nrow(t1Week), nrow(t1Week))
    # Output as list item within month:
    weekList[[j]] <- t1Week
  }
  # Bind week lists into 3 data frames that sit in the month list:
  t1List[[i]] <- bind_rows(weekList)
}

# Modify the sample events to represent sample event order within a given month rather than events within a week:

t1Frame <- bind_rows(t1List) %>%
  arrange(month, week, event) %>%
  group_by(month) %>%
  mutate(sOrder = 1:52)
  

# View the frame to double-check whether there are 4-1 records for a given site:

View(t1Frame)

# Remove unnecessary columns:

sampleOutput <- t1Frame %>%
  select(site, month, sOrder)
  


