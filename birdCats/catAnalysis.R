# =================================================================================*
# ---- SET-UP ----
# =================================================================================*
library(unmarked); library(dplyr); library(tidyr)

setwd('/Users/bsevans/Desktop/gits/birdCats/birdCats/')
list.files()

options(stringsAsFactors = F)

# ---------------------------------------------------------------------------------*
# ---- Get data ----
# ---------------------------------------------------------------------------------*

catCam <- read.csv('catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT))

catSites <- read.csv('catSites.csv') %>%
  tbl_df

catTransect <- read.csv('catDataTransect.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT))

catIncidental <- read.csv('catDataIncidental.csv') %>%
  tbl_df %>%
  filter(!is.na(CAT))

# =================================================================================*
# ---- TRANSECT DATA: DISTANCE SAMPLING ----
# =================================================================================*

# Create an unmarked frame of distance data for the transect counts

catTransectUmf <- formatDistData(catTransect, 
                                distCol="DISTANCE",
                                transectNameCol="site", 
                                dist.breaks=seq(0, 50, by = 5))
