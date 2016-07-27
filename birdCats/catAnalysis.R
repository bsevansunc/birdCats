# =================================================================================*
# ---- SET-UP ----
# =================================================================================*
library(unmarked); library(dplyr); library(tidyr)

# setwd('/Users/bsevans/Desktop/gits/birdCats/birdCats/') # Macbook -- B
# setwd('C:/Users/Brian/Desktop/gits/birdCats') # Office Windows  -- B


list.files()

options(stringsAsFactors = F)

# ---------------------------------------------------------------------------------*
# ---- Get data ----
# ---------------------------------------------------------------------------------*

catDataActivity <- read.csv('catDataActivity.csv') %>%
  tbl_df

catCam <- read.csv('catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(species))

catSites <- read.csv('catSiteData.csv') %>%
  tbl_df %>%
  select(site, imp)

catTransect <- read.csv('catDataTransect.csv') %>%
  tbl_df %>%
  filter(!is.na(COUNT)) %>%
  left_join(catDataActivity %>%
              filter(ACTIVITY == 'Transect'),
            by = c('SITE', 'VISIT')
            ) %>%
  # filter(!is.na(DISTANCE)) %>%
  mutate(DISTANCE = as.numeric(DISTANCE),
         VISIT = as.factor(VISIT)) %>%
  arrange(SITE)

catIncidental <- read.csv('catDataIncidental.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT))

# =================================================================================*
# ---- TRANSECT DATA: DISTANCE SAMPLING ----
# =================================================================================*

# Create an unmarked frame of distance data for the transect counts

catTransectUmf <- formatDistData(
  data.frame(catTransect), 
  distCol="DISTANCE",
  transectNameCol="SITE", 
  dist.breaks=seq(0,50, by = 5) # ,
  # occasionCol = 'VISIT'
  )

covs <- left_join(
  catTransect,
  read.csv('llData.csv') %>%
    tbl_df %>%
    arrange(site) %>%
    dplyr::rename(SITE = site),
  by = 'SITE'
  ) %>%
  select(-c(VISIT:TIME)) %>%
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
