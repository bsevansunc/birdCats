# =================================================================================*
# ---- SET-UP ----
# =================================================================================*
library(unmarked); library(dplyr); library(tidyr); library(camtrapR); library(ggplot2)

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
  dist.breaks=seq(0,50, by = 5),
  occasionCol = 'visit'
  )

# Get abundance covariates

sitesWithCovs <- left_join(
  catTransect,
  covs,
  by = 'site'
  ) %>%
  select(-c(visit:dew)) %>%
  distinct


# Get detection covariates

sitesWithObsCovs <- catTransect %>%
  select(-c(species:date)) %>%
  distinct


# Function that transposes covariates to wider format

transposeCovariate <- function(data, covariate){
  sites <- unique(data$site)
  transMat <- matrix(ncol = 4, nrow = length(sites))
  for(i in 1:length(sites)){
    dataSub <- data %>% filter(site == sites[i])
    transMat[i,] <- t(dataSub[covariate])
    }
  transDf <- transMat %>% data.frame
  return(transDf)
  }


# Combining wide detection covariate frames:

transTime <- transposeCovariate(sitesWithObsCovs, 'time')
transTemp <- transposeCovariate(sitesWithObsCovs, 'temp')
transDew <- transposeCovariate(sitesWithObsCovs, 'dew')

longCovs <- cbind.data.frame(transTime,transTemp, transDew)
colnames(longCovs) <- c(rep('time', 4), rep('temp', 4), rep('dew', 4))

longSitesWithObsCovs <- sitesWithObsCovs %>%
  select(site) %>%
  unique %>%
  cbind.data.frame(longCovs)


# Create unmarkedFrameGDS object for gdistsamp

umfWithCovs <- unmarkedFrameGDS(
  y =              as.matrix(catTransectUmf),
  siteCovs =       data.frame(sitesWithCovs),
  numPrimary =     4,                              
  yearlySiteCovs = longSitesWithObsCovs, #data.frame(sitesWithObsCovs),    # For use with unmarkedFrameGDS
  survey =         'line',
  dist.breaks=     seq(0,50, by = 5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
  )


# ---------------------------------------------------------------------------------*
# ----Transect model fitting----
# ---------------------------------------------------------------------------------*


gDensityNull <- gdistsamp(~1, ~1, ~1, umfWithCovs)

gDensityImp <- gdistsamp(~imp, ~1, ~1, umfWithCovs)

gDetTemp <- gdistsamp(~1, ~1, ~temp, umfWithCovs)

gDensityEduC_detDew <- gdistsamp(~eduC, ~1, ~dew, umfWithCovs)


# Null:

densityNull <- distsamp(~1 ~1, umfWithCovs)


# Global:

densityGlobal <-distsamp(~1 ~can +age + marred + eduC, umfWithCovs)


# Single covs:

densityHumanDensity <- distsamp(~1 ~hDensity, umfWithCovs) # Model did not converge
densityImp <- distsamp(~1 ~imp, umfWithCovs)
densityAge <- distsamp(~1 ~age, umfWithCovs)
densityIncome <- distsamp(~1 ~medianIncome, umfWithCovs) # NaNs produced
densityEduHS <- distsamp(~1 ~eduHS, umfWithCovs)
densityEduC <- distsamp(~1 ~eduC, umfWithCovs)
densityCan <- distsamp(~1 ~can, umfWithCovs)
densityMar <- distsamp(~1 ~marred, umfWithCovs)



# Additive:

densityHDensityAge <- distsamp(~1 ~hDensity + age, umfWithCovs) # Model did not converge
densityMarAge <- distsamp(~1 ~marred + age, umfWithCovs)
densityCanImp <- distsamp(~1 ~can + imp, umfWithCovs)
densityIncomeAge <- distsamp(~1 ~medianIncome + age, umfWithCovs) # NaNs produced
densityAgeEduC <- distsamp(~1 ~eduC + age, umfWithCovs)
densityCanAge <- distsamp(~1 ~can + age, umfWithCovs)
densityCanEduC <- distsamp (~1 ~can + eduC, umfWithCovs)
densityMarEduC <- distsamp(~1 ~marred + eduC, umfWithCovs)

# Interaction:

densityImpIntAge <- distsamp(~1 ~imp*age, umfWithCovs) # NaNs produced
densityCanIntAge <- distsamp(~1 ~can*age, umfWithCovs)
densityImpIntHDensity <- distsamp(~1 ~imp*hDensity, umfWithCovs)
densityMarIntAge <- distsamp(~1 ~marred*age, umfWithCovs)
densityCanIntHDensity <- distsamp(~1 ~can*hDensity, umfWithCovs) # Model did not converge
densityMarIntEduC <- distsamp(~1 ~marred*eduC, umfWithCovs)


# View info, example:

densityNull
logLik(densityNull)*(-2)



# Get density estimates from model

siteDensity <- predict(densityGlobal, type="state") %>%
  select(Predicted) %>%
  data.frame



# =================================================================================*
# ---- CAMERA DATA ----
# =================================================================================*

# Create detection history for each site:

umfCam <- read.csv('catCamDetection.csv') %>%
  data.frame


# Get abundance covariates for the camera-only sites

camCovs <- covs %>%
  filter(
    site != 'GERYERIMD1' &
    site != 'OLONMARDC1' &
    site != 'MISSEDDC1' &
    site != 'WOLFKARDC1' &
    site != 'WOLFAMYDC1'
  ) %>%
  data.frame


# Get detection covariates

camDetCovs <- read.csv('camDetCovs.csv') %>%
  select(site, day, tempHigh, tempLow, dewLow)


# Create and unmarkedFramePCount object for pcount

camUmfWithCovs <- unmarkedFramePCount(
  umfCam[,-1],
  siteCovs = camCovs,
  obsCovs = camDetCovs
  )


# ---------------------------------------------------------------------------------*
# ----Cam model fitting----
#----------------------------------------------------------------------------------*


#Null

camDensityNull <- pcount(~1 ~1, camUmfWithCovs, K = 50)


# Single abundance covariates

camDensityImp <- pcount(~1 ~imp, camUmfWithCovs, K = 50)
camDensityCan <- pcount(~1 ~can, camUmfWithCovs, K = 50)
camDensityhDensity <- pcount(~1 ~hDensity, camUmfWithCovs, K = 50)
camDensityAge <- pcount(~1 ~age, camUmfWithCovs, K = 50)
camDensityIncome <- pcount(~1 ~medianIncome, camUmfWithCovs, K = 50)
camDensityEduC <- pcount(~1 ~eduC, camUmfWithCovs, K = 50)
camDensityEduHS <- pcount(~1 ~eduHS, camUmfWithCovs, K = 50)
camDensityMar <- pcount(~1 ~marred, camUmfWithCovs, K = 50)



# Get density estimates from model:

camSiteDensity <- predict(camDensityNull, type = 'state') %>%
  select(Predicted) %>%
  data.frame



# ================================================================================*
# --------- PLOT ----------
# ================================================================================*

ggplot(data = siteDensity, aes(x = Predicted)) +
  geom_histogram(bins = 10)

ggplot(data = camSiteDensity, aes(x = Predicted)) +
  geom_histogram(bins = 10)
