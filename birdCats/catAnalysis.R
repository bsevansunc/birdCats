# =================================================================================*
# ---- SET-UP ----
# =================================================================================*

# Function searches packages in installed package list, add them if you don't have them, and loads the library:

smartLibrary <- function(packageVector){
  for(i in 1:length(packageVector)){
    package <- packageVector[i]
    if(!package %in% rownames(installed.packages())){
      install.packages(packageVector[i],repos="http://cran.rstudio.com/",
                       dependencies=TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}

smartLibrary(c('unmarked', 'dplyr', 'tidyr', 'camtrapR', 'ggplot2', 'AICcmodavg'))

# setwd('/Users/bsevans/Desktop/gits/birdCats/birdCats/') # Macbook -- B
# setwd('C:/Users/Brian/Desktop/gits/birdCats') # Office Windows  -- B
# setwd('C:/Users/kbenn/Documents/GitHub/birdCats/birdCats') # Laptop -- K


# list.files()

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
  transMat <- matrix(ncol = 6, nrow = length(sites))
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
colnames(longCovs) <- c(rep('time', 6), rep('temp', 6), rep('dew', 6))

longSitesWithObsCovs <- sitesWithObsCovs %>%
  select(site) %>%
  unique %>%
  cbind.data.frame(longCovs)


# Split the obsCovs into different matrices
# as suggested on the unmarked Google Group

longCovsTime <- longCovs[,1:6] %>%
  as.matrix
colnames(longCovsTime) <- NULL


longCovsTemp <- longCovs[,7:12] %>%
  as.matrix
colnames(longCovsTemp) <- NULL


longCovsDew <- longCovs[,13:18] %>%
  as.matrix
colnames(longCovsDew) <- NULL
  


# Create unmarkedFrameGDS object for gdistsamp

gUmfWithCovs <- unmarkedFrameGDS(
  y =              as.matrix(catTransectUmf),
  siteCovs =       data.frame(sitesWithCovs),
  numPrimary =     6,
  yearlySiteCovs = list(time = longCovsTime, temp = longCovsTemp, dew = longCovsDew),
  survey =         'line',
  dist.breaks=     seq(0,50, by = 5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
  )


# Create unmarkedFrameDS object for distsamp

# umfWithCovs <- unmarkedFrameDS(
#   y =              as.matrix(catTransectUmf),
#   siteCovs =       data.frame(sitesWithCovs),               
#   survey =         'line',
#   dist.breaks=     seq(0,50, by = 5),
#   tlength =        rep(200, nrow(catTransectUmf)),
#   unitsIn =        'm'
#   )


# ---------------------------------------------------------------------------------*
# ----Transect model fitting----
# ---------------------------------------------------------------------------------*

# gdistsamp

# Null model

mNull <-  gdistsamp(
  lambdaformula = ~ 1,
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density', 
  unitsOut = 'ha', mixture="P")


# Global model

mGlobal <- gdistsamp(
  lambdaformula = ~ scale(imp)*scale(can) +
    scale(imp^2) + scale(can^2) +
    scale(age)*scale(marred) + scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density', 
  unitsOut = 'ha', mixture="P")


# Abundance--impervious^2 + median income. This is the best model

mImp2Income <- gdistsamp(
  lambdaformula = ~ scale(imp)  + scale(imp^2) + 
    scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density', unitsOut = 'ha', 
  mixture="P")


# Abundance--imp

mImp <- gdistsamp(
  lambdaformula = ~ scale(imp),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Abundance--can

mCan <- gdistsamp(
  lambdaformula = ~ scale(can),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Abundance--medianIncome

mIncome <- gdistsamp(
  lambdaformula = ~ scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Abundance--imp + can

mImpCan <- gdistsamp(
  lambdaformula = ~ scale(imp) + scale(can),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Abundance--imp^2 + can^2 + medianIncome

mImp2Can2Income <- gdistsamp(
  lambdaformula = ~ scale(imp) + scale(imp^2) + 
    scale(can) + scale(can^2) + scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density', 
  unitsOut = 'ha', mixture="P")


# Abundance--impervious^2

mImp2 <- gdistsamp(
  lambdaformula = ~ scale(imp) + scale(imp^2),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density', 
  unitsOut = 'ha', mixture="P")


# Abundance--canopy^2

mCan2 <- gdistsamp(
  lambdaformula = ~ scale(can) + scale(can^2),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density',
  unitsOut = 'ha', mixture="P")


# Abundance--canopy^2 + medianIncome

mCan2Income <- gdistsamp(
  lambdaformula = ~ scale(can) + scale(can^2) +
    scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density',
  unitsOut = 'ha', mixture="P")


# Abundance--can + medianIncome

mCanIncome <- gdistsamp(
  lambdaformula = ~ scale(can) +
    scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density',
  unitsOut = 'ha', mixture="NB")


# Abundance--imp + medianIncome

mImpIncome <- gdistsamp(
  lambdaformula = ~ scale(imp) +
    scale(medianIncome),
  phiformula = ~ scale(time),
  pformula = ~ scale(dew) * scale(temp) + scale(time),
  data = gUmfWithCovs,
  keyfun = "halfnorm", output = 'density',
  unitsOut = 'ha', mixture="NB")


# -----------------*
# ---- Density ----
# -----------------*

transDensity <- predict(mImp2Income, type="lambda", appenddata = TRUE)


# ----------------------------------------------------------------*
# ------ Predict based on new data ------
# ----------------------------------------------------------------*

nd <- data.frame(imp = 0:100, 
                 medianIncome = 106227,
                 time = 533.5,
                 dew = 68,
                 temp = 79)

predictImp <- predict(mImp2Income, newdata = nd, type = 'lambda', appenddata = TRUE)


nd2 <- data.frame(imp = median(as.vector(covs$imp)), 
                 medianIncome = seq(49800, 240600, by = 100),
                 time = median(as.vector(sitesWithObsCovs$time)),
                 dew = median(as.vector(sitesWithObsCovs$dew)),
                 temp = median(as.vector(sitesWithObsCovs$temp)))

predictInc <- predict(mImp2Income, newdata = nd2, type = 'lambda', appenddata = TRUE)



# Get AICc values from models:

AICc(mNull, return.K = FALSE)
AICc(mGlobal, return.K = FALSE)
AICc(mCan, return.K = FALSE)
AICc(mImp, return.K = FALSE)
AICc(mIncome, return.K = FALSE)
AICc(mImpCan, return.K = FALSE)
AICc(mImp2, return.K = FALSE) # Second-best
AICc(mCan2, return.K = FALSE)
AICc(mCanIncome, return.K = FALSE)
AICc(mImpIncome, return.K = FALSE)
AICc(mImp2Income, return.K = FALSE) # This is the best model
AICc(mCan2Income, return.K = FALSE)
AICc(mImp2Can2Income, return.K = FALSE)


# Get log-likelihoods:

extractLL(mNull)*(-2)
extractLL(mCan)*(-2)
extractLL(mImp)*(-2)
extractLL(mImpCan)*(-2)
extractLL(mCan2)*(-2)
extractLL(mImp2)*(-2)
extractLL(mCan2Income)*(-2)
extractLL(mImp2Income)*(-2)
extractLL(mImp2Can2Income)*(-2)
extractLL(mIncome)*(-2)
extractLL(mCanIncome)*(-2)
extractLL(mImpIncome)*(-2)



# ----------------------------------------------------------------*
# ------ Visualize variables ------
# ----------------------------------------------------------------*

transSites <- covs %>%
  select(site, imp, medianIncome)

transSiteDensity <- cbind.data.frame(transSites, transDensity$Predicted)
imp <- transSiteDensity$imp
inc <- transSiteDensity$medianIncome
dens <- transSiteDensity$Predicted

plot(imp, dens)
plot(inc, dens)




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

camNull <- pcount(
  formula = ~scale(dewLow) ~1, 
  data= camUmfWithCovs,
  K = 50)



# Global

camImp2Can2Income <- pcount(
  formula = ~scale(dewLow) 
    ~scale(can) + scale(can^2) + scale(imp) + scale(imp^2) + scale(medianIncome), 
  data= camUmfWithCovs,
  K = 50)


# Single abundance covariates

camImp <- pcount(
  formula = ~dewLow ~imp, 
  data = camUmfWithCovs,
  K = 50)

camCan <- pcount(
  formula = ~dewLow ~can, 
  data = camUmfWithCovs, 
  K = 50)

camIncome <- pcount(
  formula = ~scale(dewLow) ~scale(medianIncome),
  data = camUmfWithCovs, 
  K = 50)



# Additive

camCanIncome <- pcount(
  formula = ~scale(dewLow) ~scale(can) + scale(medianIncome),
  data = camUmfWithCovs,
  K = 50)

camImpIncome <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(medianIncome),
  data = camUmfWithCovs,
  K = 50)

camImpCan <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(can),
  data = camUmfWithCovs,
  K = 50)


# Quadratic

camImp2 <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(imp^2),
  data = camUmfWithCovs,
  K = 50)

# This is the best model: 
camCan2 <- pcount(
  formula = ~scale(dewLow) ~scale(can) + scale(can^2),
  data = camUmfWithCovs,
  K = 50)

camImp2Income <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(imp^2) + scale(medianIncome),
  data = camUmfWithCovs,
  K = 50)

camCan2Income <- pcount(
  formula = ~scale(dewLow) ~scale(can) + scale(can^2) + scale(medianIncome),
  data = camUmfWithCovs,
  K = 50)

camCan2Imp <- pcount(
  formula = ~scale(dewLow) ~scale(can) + scale(can^2) + scale(imp),
  data = camUmfWithCovs,
  K = 50)

camImp2Can <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(imp^2) + scale(can),
  data = camUmfWithCovs,
  K = 50)

camImp2Can2 <- pcount(
  formula = ~scale(dewLow) ~scale(imp) + scale(imp^2) + scale(can) + scale(can^2),
  data = camUmfWithCovs,
  K = 50)


# -----------------*
# ---- Density ----
# -----------------*

camDensity <- predict(camCan2, type = 'state') %>%
  select(Predicted) %>%
  data.frame


# -------------------------------------------------*
# ------ Visualize variables ------
# -------------------------------------------------*

camSites <- camCovs %>%
  select(site, can)

camSiteDensity <- cbind.data.frame(camSites, camDensity)
can <- camSiteDensity$can
camDens <- camSiteDensity$Predicted

plot(can, camDens)



# Get AICc scores

AICc(camNull, return.K = FALSE)
AICc(camCan, return.K = FALSE)
AICc(camImp, return.K = FALSE)
AICc(camIncome, return.K = FALSE)
AICc(camCan2, return.K = FALSE)
AICc(camImp2, return.K = FALSE)
AICc(camImpCan, return.K = FALSE)
AICc(camImpIncome, return.K = FALSE)
AICc(camCanIncome, return.K = FALSE)
AICc(camCan2Income, return.K = FALSE)
AICc(camImp2Income, return.K = FALSE)
AICc(camImp2Can2Income, return.K = FALSE)
AICc(camImp2Can2, return.K = FALSE)
AICc(camImp2Can, return.K = FALSE)
AICc(camCan2Imp, return.K = FALSE)


# Get log-likelihood:

logLik(camNull)*(-2)
logLik(camCan)*(-2)
logLik(camImp)*(-2)
logLik(camIncome)*(-2)
logLik(camCan2)*(-2)
logLik(camImp2)*(-2)
logLik(camImpCan)*(-2)
logLik(camImpIncome)*(-2)
logLik(camCanIncome)*(-2)
logLik(camCan2Income)*(-2)
logLik(camImp2Income)*(-2)
logLik(camImp2Can2)*(-2)
logLik(camImp2Can)*(-2)
logLik(camCan2Imp)*(-2)
logLik(camImp2Can2Income)*(-2)





# ================================================================================*
# --------- PLOT ----------
# ================================================================================*

# THIS WILL NOT WORK TEMPORARILY--I have to combine the trans and cam
# density estimates better before this will work again. Density estimates 
# are plotted against individual variables at the end of the trans and cam
# sections above.


# Combine density estimates:

# transSites <- catTransect %>%
#   select(site) %>%
#   unique
# 
# camSites <- umfCam %>%
#   select(site) %>%
#   unique
# 
# transSiteDensity <- cbind.data.frame(transSites, transDensity)
# length(camDensity) <- 53
# camDensity <- camDensity %>% data.frame
# 
# 
# siteDensityUntidy <- cbind.data.frame(transSiteDensity, camSiteDensity)
# 
# colnames(siteDensityUntidy) <- c('site', 'trans', 'cam')
# 
# # Tidy up the data:
# 
# siteDensity <- gather(data = siteDensityUntidy,
#                       key = method,
#                       value = density,
#                       trans:cam) %>%
#   na.omit %>%
#   group_by(method)
# 
# densitySumm <- siteDensity %>%
#   summarise(sample_size = length(method), 
#             mean = mean(density), 
#             sd = sd(density),
#             se = sd(density)/sqrt(length(method)))
# 
# 
# densPlot <- ggplot(data = densitySumm, aes(x = method, y = mean))+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0, size = 1)+
#   geom_point(fill = 'red', color = 'black', shape = 21, size = 5)+
#   scale_x_discrete('Method', labels = c('Camera', 'Transect'))+
#   scale_y_continuous('Mean density (cats/ha)', limits = c(0, 2.2))+
#   theme(panel.grid = element_blank(), 
#         axis.line.x = element_line(linetype='solid', color='black'),
#         axis.line.y = element_line(linetype='solid', color='black'))
