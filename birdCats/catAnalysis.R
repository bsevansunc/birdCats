# =================================================================================*
# ---- SET-UP ----
# =================================================================================*
library(unmarked); library(dplyr); library(tidyr)

# setwd('/Users/bsevans/Desktop/gits/birdCats/birdCats/') # Macbook -- B
# setwd('C:/Users/Brian/Desktop/gits/birdCats') # Office Windows  -- B
# setwd('C:/Users/Kevin/Documents/GitHub/birdCats/birdCats') # Laptop -- K


list.files()

options(stringsAsFactors = F)

# ---------------------------------------------------------------------------------*
# ---- Get data ----
# ---------------------------------------------------------------------------------*

catSiteActivity <- read.csv('catDataActivity.csv') %>%
  tbl_df

catCam <- read.csv('catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(species))

catSites <- read.csv('catSiteData.csv') %>%
  tbl_df %>%
  select(site, imp)

catTransect <- read.csv('catDataTransect.csv') %>%
  tbl_df %>%
  filter(!is.na(count)) %>%
  filter(species == 'cat') %>%
  left_join(catSiteActivity %>%
              filter(activity == 'transect'),
            by = c('site', 'visit')
            ) %>%
  # filter(!is.na(distance)) %>%
  mutate(distance = as.numeric(distance),
         visit = as.factor(visit)) %>%
  arrange(site)

catIncidental <- read.csv('catDataIncidental.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT))

# =================================================================================*
# ---- TRANSECT DATA: DISTANCE SAMPLING ----
# =================================================================================*

# Create an unmarked frame of distance data for the transect counts

catTransectUmf <- formatDistData(
  data.frame(catTransect), 
  distCol="distance",
  transectNameCol="site", 
  dist.breaks=seq(0,50, by = 5) # ,
  # occasionCol = 'visit'
  )

covs <- left_join(
  catTransect,
  read.csv('llData.csv') %>%
    tbl_df %>%
    arrange(site),
  by = 'site'
  ) %>%
  select(-c(visit:time)) %>%
  distinct
  
umfWithCovs <- unmarkedFrameDS(
  y = as.matrix(catTransectUmf),
  siteCovs = data.frame(covs),
  survey = 'line',
  dist.breaks=seq(0,50, by = 5),
  tlength = rep(200, nrow(catTransectUmf)),
  unitsIn = 'm'
  )

umfWithCovs %>% summary

# 
