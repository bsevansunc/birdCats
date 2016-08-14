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
  mutate(medianIncomeAdj = medianIncome/1000) %>%
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
  dist.breaks=seq(0,50, by = 5)#,
#  occasionCol = 'visit'
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

# umfWithCovs <- unmarkedFrameGDS(
#   y =              as.matrix(catTransectUmf),
#   siteCovs =       data.frame(sitesWithCovs),
#   numPrimary =     1,                              
#   yearlySiteCovs = data.frame(longSitesWithObsCovs[,-1]),
#   survey =         'line',
#   dist.breaks=     seq(0,50, by = 5),
#   tlength =        rep(200, nrow(catTransectUmf)),
#   unitsIn =        'm'
#   )


# Create unmarkedFrameDS object for distsamp

umfWithCovs <- unmarkedFrameDS(
  y =              as.matrix(catTransectUmf),
  siteCovs =       data.frame(sitesWithCovs),               
  survey =         'line',
  dist.breaks=     seq(0,50, by = 5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
  )


# ---------------------------------------------------------------------------------*
# ----Transect model fitting----
# ---------------------------------------------------------------------------------*

# Not working

# gDensityNull <- gdistsamp(~1, ~1, ~1, umfWithCovs)
# 
# gDensityImp <- gdistsamp(~imp, ~1, ~1, umfWithCovs)
# 
# gDetTemp <- gdistsamp(~1, ~1, ~temp, umfWithCovs)
# 
# gDensityEduC_detDew <- gdistsamp(~eduC, ~1, ~dew, umfWithCovs)
# 
# gDetDewTemp <- gdistsamp(~1, ~1, ~dew+temp, umfWithCovs)
# 
# gDensityImp_detDewTemp <- gdistsamp(~imp, ~1, ~dew+temp, umfWithCovs)




# Null:

densityNull <- distsamp(~1 ~1, umfWithCovs)


# Global:

densityGlobal <-distsamp(
  ~1 ~can + hDensity + age + marred + eduC + medianIncomeAdj, umfWithCovs
  )


# Single covs:

densityHDensity <- distsamp(
  ~1 ~hDensity, umfWithCovs, control=list(
    maxit=1000,trace=TRUE, REPORT=1
    )
  )

densityImp <- distsamp(~1 ~imp, umfWithCovs)

densityAge <- distsamp(~1 ~age, umfWithCovs)

densityIncome <- distsamp(~1 ~medianIncomeAdj, umfWithCovs) # NaNs produced

densityEduHS <- distsamp(~1 ~eduHS, umfWithCovs)

densityEduC <- distsamp(~1 ~eduC, umfWithCovs)

densityCan <- distsamp(~1 ~can, umfWithCovs)

densityMar <- distsamp(~1 ~marred, umfWithCovs)



# Additive:

densityHDensityAge <- distsamp(
  ~1 ~hDensity + age, umfWithCovs, control=list(
    maxit=1000,trace=TRUE, REPORT=1
    )
  )

densityHDensityEduC <- distsamp(
  ~1 ~hDensity + eduC, umfWithCovs, control=list(
    maxit=1700,trace=TRUE, REPORT=1
    )
  )

densityHDensityCan <- distsamp(
  ~1 ~hDensity + can, umfWithCovs, control=list(
    maxit=1700,trace=TRUE, REPORT=1
  )
)

densityHDensityCanEduC <- distsamp(
  ~1 ~hDensity + can + eduC, umfWithCovs, control=list(
    maxit=1700,trace=TRUE, REPORT=1
  )
)

densityMarAge <- distsamp(~1 ~marred + age, umfWithCovs)

densityCanImp <- distsamp(~1 ~can + imp, umfWithCovs)

densityIncomeAge <- distsamp(~1 ~medianIncomeAdj + age, umfWithCovs)

densityAgeEduC <- distsamp(~1 ~eduC + age, umfWithCovs)

densityCanAge <- distsamp(~1 ~can + age, umfWithCovs)

densityCanEduC <- distsamp (~1 ~can + eduC, umfWithCovs)

densityMarEduC <- distsamp(~1 ~marred + eduC, umfWithCovs)

densityCanIncome <- distsamp(~1 ~can + medianIncomeAdj, umfWithCovs)

densityIncomeEduC <- distsamp(~1 ~medianIncomeAdj + eduC, umfWithCovs)

densityIncomeEduCCan <- distsamp(~1 ~medianIncomeAdj + eduC + can, umfWithCovs)




# Interaction:

densityImpIntAge <- distsamp(~1 ~imp*age, umfWithCovs) # NaNs produced

densityCanIntAge <- distsamp(~1 ~can*age, umfWithCovs)

densityImpIntHDensity <- distsamp(~1 ~imp*hDensity, umfWithCovs)

densityMarIntAge <- distsamp(~1 ~marred*age, umfWithCovs)

densityCanIntHDensity <- distsamp(
  ~1 ~can*hDensity, umfWithCovs, control=list(
    maxit=700,trace=TRUE, REPORT=1
    )
  ) 

densityMarIntEduC <- distsamp(~1 ~marred*eduC, umfWithCovs)


# View info, example:

densityNull
logLik(densityNull)*(-2)



# Get density estimates from model

transSiteDensity <- predict(densityGlobal, type="state") %>%
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
  siteCovs = camCovs#,
#  obsCovs = camDetCovs
  )


# ---------------------------------------------------------------------------------*
# ----Cam model fitting----
#----------------------------------------------------------------------------------*


#Null

camDensityNull <- pcount(~1 ~1, camUmfWithCovs, K = 50)

camDetDew <- pcount(~dewLow ~1, camUmfWithCovs, K = 50)
camDetTempHigh <- pcount(~tempHigh ~1, camUmfWithCovs, K = 50)


# Global

camDensityGlobal <- pcount(
  ~1 ~can+hDensity+medianIncomeAdj+eduC+marred+age, camUmfWithCovs, K = 50
  )


# Single abundance covariates

camDensityImp <- pcount(~1 ~imp, camUmfWithCovs, K = 50)
camDensityCan <- pcount(~1 ~can, camUmfWithCovs, K = 50)
camDensityhDensity <- pcount(~1 ~hDensity, camUmfWithCovs, K = 50)
camDensityAge <- pcount(~1 ~age, camUmfWithCovs, K = 50)
camDensityIncome <- pcount(~1 ~medianIncomeAdj, camUmfWithCovs, K = 50)
camDensityEduC <- pcount(~1 ~eduC, camUmfWithCovs, K = 50)
camDensityEduHS <- pcount(~1 ~eduHS, camUmfWithCovs, K = 50)
camDensityMar <- pcount(~1 ~marred, camUmfWithCovs, K = 50)


# Additive

camDensityCanEduC <- pcount(~1 ~can + eduC, camUmfWithCovs, K = 50)
camDensityCanIncome <- pcount(~1 ~can + medianIncomeAdj, camUmfWithCovs, K = 50)
camDensityAgeEduC <- pcount(~1 ~age + eduC, camUmfWithCovs, K = 50)
camDensityAgeCan <- pcount(~1 ~age + can, camUmfWithCovs, K = 50)

camDensityAgeCanEduC <- pcount(~1 ~age + can + eduC, camUmfWithCovs, K = 50)





# Get density estimates from model:

camSiteDensity <- predict(camDensityGlobal, type = 'state') %>%
  select(Predicted) %>%
  as.matrix



# ================================================================================*
# --------- PLOT ----------
# ================================================================================*

# Combine density estimates:

transSites <- catTransect %>%
  select(site) %>%
  unique

# camSites <- umfCam %>%
#   select(site) %>%
#   unique

transSiteDensity <- cbind.data.frame(transSites, transSiteDensity)
length(camSiteDensity) <- 53
camSiteDensity <- camSiteDensity %>% data.frame


siteDensityUntidy <- cbind.data.frame(transSiteDensity, camSiteDensity)

colnames(siteDensityUntidy) <- c('site', 'trans', 'cam')

# Tidy up the data:

siteDensity <- gather(data = siteDensityUntidy,
                      key = method,
                      value = density,
                      trans:cam) %>%
  na.omit %>%
  group_by(method)

densitySumm <- siteDensity %>%
  summarise(sample_size = length(method), 
            mean = mean(density), 
            sd = sd(density),
            se = sd(density)/sqrt(length(method)))


ggplot(data = densitySumm, aes(x = method, y = mean))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0, size = 1)+
  geom_point(fill = 'red', color = 'black', shape = 21, size = 4)+
  scale_x_discrete('Method', labels = c('Camera', 'Transect'))+
  scale_y_continuous('Mean density (per hectare)', limits = c(0, 1.2))+
  theme(panel.grid = element_blank(), 
        axis.line.x = element_line(linetype='solid', color='black'),
        axis.line.y = element_line(linetype='solid', color='black'))
