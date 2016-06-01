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
  mutate(month = rep(1:3, length(unique(site)))) %>%
  group_by(site, month) %>%
  mutate(week = sample(1:4, 1)) %>%
  ungroup %>%
  arrange(site) 

tList <- vector('list', length(catSites))

for(i in 1:length(tList)){
  siteItem <- filter(t1, site == catSites[[i]]) %>%
    arrange(month)
  for(j in 2:3){
    if(siteItem[j,3] == 1 & siteItem[j-1,3] == 4){
      siteItem[j, 3] <- sample(2:4, 1)
    }
  }
  tList[[i]] <- siteItem
}

t1 <- bind_rows(tList)

# Make a data frame of sample events by week

t1List <- vector('list', length = 3)

for(i in 1:3){
  t1Month <- t1 %>%
    filter(month == i)
  weekList <- vector('list', length = length(unique(t1Month$week)))
  for(j in 1:length(weekList)){
    t1Week <- filter(t1Month, week == j)
    t1Week$event <- sample(1:nrow(t1Week), nrow(t1Week))
    weekList[[j]] <- t1Week
  }
  t1List[[i]] <- bind_rows(weekList)
}

# modify the sample events to represent sample event order within a given month rather than events within a week:

t1Frame <- bind_rows(t1List) %>%
  arrange(month, week, event) %>%
  group_by(month) %>%
  mutate(sOrder = 1:52)
  select(site, month, sOrder)

# View the frame:

View(t1Frame)


