# =================================================================================*
# ---- SET-UP ----
# =================================================================================*
library(unmarked); library(dplyr); library(tidyr); library(camtrapR)

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
  filter(!is.na(species)) %>%
  filter(!is.na(cameraID))

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



# =================================================================================*
# ---- CAMERA DATA ----
# =================================================================================*



# Playing around to see how recordTable() and detectionHistory() format things:


# Creating a dataframe camtrapR will like for use in cameraOperation():
#     Problem--need setup and retrieval dedicated cols


cameraActivity <- catSiteActivity %>%
                    filter(activity == "setup" | activity == "takedown")


# Did not work:

# camSetup <- cameraActivity %>%
#               filter(activity == 'setup') %>%
#               select(site, date)

# camTakedown <- cameraActivity %>%
#                 filter(activity == 'takedown') %>%
#                 select(site, date)

# camReady <- left_join(camSetup, 
#                       camTakedown,
#                       by = 'site'
#                       )
 


data(recordTableSample)

view(recordTableSample)


# Not ready yet:

# cameraOperation(cameraActivity,
#                 stationCol = "site",
#                 setupCol = "")


# Not ready yet:

# detectionHistory(recordTableSample,
#                  species = "PBE",
#                  )


