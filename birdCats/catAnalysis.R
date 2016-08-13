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
  filter(visit <= 4) %>%
  filter(!is.na(count)) %>%
  filter(species == 'cat') %>%
  left_join(catSiteActivity %>%
              filter(activity == 'transect'),
            by = c('site', 'visit')
            ) %>%
  mutate(distance = as.numeric(distance),
         visit = as.factor(visit)) %>%
  arrange(site)

covs <- read.csv('covariateData.csv') %>%
  tbl_df %>%
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

sitesWithCovs <- left_join(
  catTransect,
  covs,
  by = 'site'
  ) %>%
  select(-c(visit:dew)) %>%
  distinct

sitesWithObsCovs <- catTransect %>%
  select(-c(species:date)) %>%
  distinct

# There is a problem here:

umfWithCovs <- unmarkedFrameGDS(
  y =              as.matrix(catTransectUmf),
  siteCovs =       data.frame(sitesWithCovs),
  numPrimary =     4,
  yearlySiteCovs = data.frame(sitesWithObsCovs),
  survey =         'line',
  dist.breaks=     seq(0,50, by = 5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
  )

summary(umfWithCovs)

# Fit models:

# This does not work:

gDensityNull <- gdistsamp(~1, ~1, ~1, umfWithCovs)



# Null:

densityNull <- distsamp(~1 ~1, umfWithCovs)


# Single covs:

densityHumanDensity <- distsamp(~1 ~hDensity, umfWithCovs) # Model did not converge
densityImp <- distsamp(~1 ~imp, umfWithCovs)
densityAge <- distsamp(~1 ~age, umfWithCovs)
densityIncome <- distsamp(~1 ~medianIncome, umfWithCovs) # NaNs produced
densityEduHS <- distsamp(~1 ~eduHS, umfWithCovs)
densityCan <- distsamp(~1 ~can, umfWithCovs)
densityCanImp <- distsamp(~1 ~can + imp, umfWithCovs)
densityMar <- distsamp(~1 ~marred, umfWithCovs)

densityAge_detImp <- distsamp(~imp ~age, umfWithCovs)
densityNull_detImp <- distsamp(~imp ~1, umfWithCovs)
densityEduHS_detImp <- distsamp(~imp ~eduHS, umfWithCovs)



# Additive:

densityHDensityAge <- distsamp(~1 ~hDensity + age, umfWithCovs) # Model did not converge
densityMarAge <- distsamp(~1 ~marred + age, umfWithCovs)


# Interaction:

densityImpIntAge <- distsamp(~1 ~imp*age, umfWithCovs) # NaNs produced
densityImpIntHDensity <- distsamp(~1 ~imp*hDensity, umfWithCovs)
densityMarIntAge <- distsamp(~1 ~marred*age, umfWithCovs)
densityCanIntHDensity <- distsamp(~1 ~can*hDensity, umfWithCovs) # Model did not converge


# View info, example:

densityNull
logLik(densityNull)*(-2)



# =================================================================================*
# ---- CAMERA DATA ----
# =================================================================================*

# Create detection history for each site:

umfCam <- read.csv('catCamDetection.csv')

camCovs <- covs %>%
  filter(
    site != 'GERYERIMD1' &
    site != 'OLONMARDC1' &
    site != 'MISSEDDC1' &
    site != 'WOLFKARDC1' &
    site != 'WOLFAMYDC1'
  )






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


# Package example:

data(recordTableSample)


# define image directory
wd_images_ID <- system.file("pictures/sample_images", package = "camtrapR")
# load station information
data(camtraps)
# create camera operation matrix
camop_no_problem <- cameraOperation(CTtable      = camtraps,
                                    stationCol   = "Station",
                                    setupCol     = "Setup_date",
                                    retrievalCol = "Retrieval_date",
                                    hasProblems  = FALSE,
                                    dateFormat   = "%d/%m/%Y"
)

# compute detection history for a species

# with effort
DetHist2 <- detectionHistory(recordTable          = recordTableSample,
                             camOp                = camop_no_problem,
                             stationCol           = "Station",
                             speciesCol           = "Species",
                             recordDateTimeCol    = "DateTimeOriginal",
                             species              = "VTA",
                             occasionLength       = 7,
                             day1                 = "station",
                             datesAsOccasionNames = FALSE,
                             includeEffort        = TRUE,
                             scaleEffort          = FALSE,
                             timeZone             = "Asia/Kuala_Lumpur"
)


DetHist2[[1]]  # detection history





