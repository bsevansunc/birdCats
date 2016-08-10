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
  select(site)

catTransect <- read.csv('catDataTransect.csv') %>%
  tbl_df %>%
  mutate(site = ifelse(site == 'OBRICHMD1', 'OBRICHRMD1', site)) %>%
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

covs <- read.csv('covariateData.csv') %>%
  tbl_df %>%
  mutate(site = ifelse(site == 'OBRICHMD1', 'OBRICHRMD1', site)) %>%
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

covs_df <- left_join(
  catTransect,
  covs,
  by = 'site'
  ) %>%
  select(-c(visit:time)) %>%
  distinct
  
umfWithCovs <- unmarkedFrameDS(
  y = as.matrix(catTransectUmf),
  siteCovs = data.frame(covs_df),
  survey = 'line',
  dist.breaks=seq(0,50, by = 5),
  tlength = rep(200, nrow(catTransectUmf)),
  unitsIn = 'm'
  )

umfWithCovs %>% summary


# Fit models:

densityNull <- distsamp(~1~1,umfWithCovs)
densityHumanDensity <- distsamp(~1~hDensity, umfWithCovs)
densityImp <- distsamp(~1~imp, umfWithCovs)
densityHumanDensityAge <- distsamp(~1~hDensity+age, umfWithCovs)
densityAge <- distsamp(~1~age, umfWithCovs)
densityImpIntAge <- distsamp(~1~imp*age, umfWithCovs)
densityImpIntHDensity <- distsamp(~1~imp*hDensity, umfWithCovs)
densityIncome <- distsamp(~1~medianIncome, umfWithCovs)
densityEdu <- distsamp(~1 ~eduHS, umfWithCovs)



# =================================================================================*
# ---- CAMERA DATA ----
# =================================================================================*



# Playing around to see how recordTable() and detectionHistory() format things:




# Get ready to create camera operation matrix:


camOperation <- catSiteActivity %>%
              filter(activity == 'setup') %>%
              select(site, date)

camTakedown <- catSiteActivity %>%
                filter(activity == 'takedown') %>%
                select(site, date)

camOperation[,3] <- camTakedown[,2]


colnames(camOperation) <- c('site', 'setup', 'takedown')


camID <- read.csv('camID.csv')


camOperation <- camOperation %>%
  cbind(camID)



# Create camera operation matrix:

camOpMatrix <- cameraOperation(camOperation,
                stationCol = 'site',
                cameraCol = 'cameraID',
                setupCol = 'setup',
                retrievalCol = 'takedown',
                byCamera = FALSE,
                allCamsOn = FALSE,
                camerasIndependent = TRUE
                )


# Not ready yet:

# detectionHistory(recordTableSample,
#                  species = "PBE",
#                  camOp = camOpMatrix,
#                  stationCol = 'Station',
#                  speciesCol = 'Species',
#                  recordDateTimeCol = "DateTimeOriginal",
#                  occasionLength = 1,
#                  day1 = 'station',
#                  timeZone = "America/New_York",
#                  includeEffort = FALSE
#                  )


# Sample:

data(recordTableSample)

