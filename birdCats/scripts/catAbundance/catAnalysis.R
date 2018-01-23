# =================================================================================*
# --------------------------------- Set-up ----------------------------------------
# =================================================================================*

# Function searches packages in installed package list,
# add them if you don't have them, and loads the library:

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

smartLibrary(
  c('unmarked', 'dplyr', 'tidyr', 'camtrapR', 'ggplot2', 'AICcmodavg','MuMIn')
  )


options(stringsAsFactors=F)

# ---------------------------------------------------------------------------------*
# ---- Get data ----
# ---------------------------------------------------------------------------------*

catSiteActivity <- read.csv('data/catDataActivity.csv') %>%
  tbl_df

catCam <- read.csv('data/catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(species)) %>%
  filter(!is.na(cameraID))

catSites <- read.csv('data/catSiteData.csv') %>%
  tbl_df %>%
  select(site)

catTransect <- read.csv('data/catDataTransect.csv') %>%
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

eduData <- read.csv('data/eduHS.csv')

covs <- read.csv('data/covariateData.csv') %>%
  select(-c(eduHS, eduC)) %>%
  left_join(eduData) %>%
  tbl_df %>%
  arrange(site)


# =================================================================================*
# ---------------------- Set up transect sampling data ----------------------------
# =================================================================================*

# Create an unmarked frame of distance data for the transect counts

catTransectUmf <- formatDistData(
  data.frame(catTransect), 
  distCol="distance",
  transectNameCol="site", 
  dist.breaks=seq(0, 50, by=5),
  occasionCol='visit'
  )

# Get abundance covariates

sitesWithCovs <- left_join(catTransect, covs, by='site') %>%
  select(-c(visit:dew)) %>%
  distinct


# Get detection covariates

sitesWithObsCovs <- catTransect %>%
  select(-c(species:date)) %>%
  distinct


# Function that transposes covariates to wider format

transposeCovariate <- function(data, covariate){
  sites <- unique(data$site)
  transMat <- matrix(ncol=6, nrow=length(sites))
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
  yearlySiteCovs = list(time=longCovsTime, temp=longCovsTemp, dew=longCovsDew),
  survey =         'line',
  dist.breaks=     seq(0, 50, by=5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
  )

# Square and scale variables

gUmfWithCovs@siteCovs <- gUmfWithCovs@siteCovs %>%
  mutate(
    can2 = can^2,
    can = scale(can)[,1],
    can2 = scale(can2)[,1],
    imp2 = imp^2,
    imp = scale(imp)[,1],
    imp2 = scale(imp2)[,1],
    medianIncome2 = medianIncome^2,
    medianIncome = scale(medianIncome)[,1],
    medianIncome2 = scale(medianIncome2)[,1],
    hDensity2 = hDensity^2,
    hDensity = scale(hDensity)[,1],
    hDensity2 = scale(hDensity2)[,1],
    age2 = age^2,
    age = scale(age)[,1],
    age2 = scale(age2)[,1],
    eduHS = scale(eduHS)[,1],
    marred = scale(marred)[,1])

gUmfWithCovs@yearlySiteCovs <- gUmfWithCovs@yearlySiteCovs %>%
  mutate(
    time2 = time^2,
    time = scale(time)[,1],
    time2 = scale(time2)[,1],
    temp = scale(temp)[,1],
    dew = scale(dew)[,1])


# ---------------------------------------------------------------------------------*
# ---- Transect null model fitting ----
# ---------------------------------------------------------------------------------*

# Rather than building models with every possible combination of phi and p
# formulas, we compare null models, find the best one, use that model's phi and
# p formulas in the abundance models.


# Null models

null1.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null2.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~dew*temp,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null3.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null4.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null5.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~temp,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null6.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~1,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

#This model does not run properly
null7.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~1,
  pformula = ~dew*temp,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null8.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~1,
  pformula = ~temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null9.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~1,
  pformula = ~time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

null10.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~1,
  pformula = ~temp,
  data = gUmfWithCovs,
  keyfun = "halfnorm",
  mixture = "NB"
  )

nullmodelList <- list(
  null1.model, null2.model, null3.model, null4.model, null5.model, null6.model,
  # null7.model,
  null8.model, null9.model, null10.model
  )

nulltnames <- c(
  'null1','null2','null3','null4','null5','null6',
  # 'null7',
  'null8','null9','null10'
  )


# Assemble AICc table

nulltable <- aictab(cand.set=nullmodelList,modnames=nulltnames)


# Examine the null models

nulltable


# =================================================================================*
# ------------------------- Set up camera sampling data ---------------------------
# =================================================================================*

# Create detection history for each site:

umfCam <- read.csv('data/catCamDetection.csv') %>%
  data.frame


# Single out transect-only sites

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1',
                 'GERYERIMD1', 'MISSEDDC1')

# Square and scale variables

camCovs <- covs %>%
  filter(!site %in% removeSites) %>%
  data.frame %>%
  mutate(
    imp2 = imp^2,
    imp = scale(imp)[,1],
    imp2 = scale(imp2)[,1],
    can2 = can^2,
    can = scale(can)[,1],
    can2 = scale(can2)[,1],
    medianIncome2 = medianIncome^2,
    medianIncome = scale(medianIncome)[,1],
    medianIncome2 = scale(medianIncome2)[,1],
    hDensity2 = hDensity^2,
    hDensity = scale(hDensity)[,1],
    hDensity2 = scale(hDensity2)[,1],
    age2 = age^2,
    age = scale(age)[,1],
    age2 = scale(age2)[,1],
    eduHS = scale(eduHS)[,1],
    marred = scale(marred)[,1])


# Get scaled detection covariates

camDetCovs <- read.csv('data/camDetCovs.csv') %>%
  filter(!site %in% removeSites) %>%
  select(site, day, tempHigh, tempLow, dewLow) %>%
  mutate(
    tempHigh = scale(tempHigh)[,1],
    tempLow = scale(tempLow)[,1],
    dewLow = scale(dewLow)[,1]
    )


# Create an unmarkedFramePCount object for pcount

camUmfWithCovs <- unmarkedFramePCount(
  umfCam[,-1],
  siteCovs=camCovs,
  obsCovs=camDetCovs)


# ---------------------------------------------------------------------------------*
# ---- Cam null model fitting ----
#----------------------------------------------------------------------------------*

# Same as transect models, we first compare nulls, then use the best model's phi
# and p formulas in future abundance models

#Null models

null1.cmodel <- pcount(
  formula = ~dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null2.cmodel <- pcount(
  formula = ~tempHigh ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null3.cmodel <- pcount(
  formula = ~tempLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null4.cmodel <- pcount(
  formula = ~tempLow+dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null5.cmodel <- pcount(
  formula = ~tempHigh+dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null6.cmodel <- pcount(
  formula = ~tempHigh*dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

null7.cmodel <- pcount(
  formula = ~tempLow*dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB',
  K = 50
  )

nullclist <- list(null1.cmodel,null2.cmodel,null3.cmodel,
                  null4.cmodel,null5.cmodel,null6.cmodel,null7.cmodel)

nullcnames <- c('null1','null2','null3','null4','null5','null6','null7')

nullctable <- aictab(nullclist,nullcnames)


# Examine the camera null models

nullctable

# =================================================================================*
# ---------------------- Model fitting and AIC table creation ---------------------
# =================================================================================*


# List of model formulas

formulalist <- list(
  ~imp, ~imp+imp2, ~can, ~can+can2, ~hDensity, ~hDensity+hDensity2, ~imp+age,
  ~imp+marred, ~imp+medianIncome, ~imp+eduHS, ~imp+imp2+age, ~imp+imp2+marred,
  ~imp+imp2+medianIncome, ~imp+imp2+eduHS, ~can+age, ~can+marred, ~can+medianIncome,
  ~can+eduHS, ~can+can2+age, ~can+can2+marred, ~can+can2+medianIncome,
  ~can+can2+eduHS, ~hDensity+age, ~hDensity+marred, ~hDensity+medianIncome,
  ~hDensity+eduHS, ~hDensity+hDensity2+age, ~hDensity+hDensity2+marred,
  ~hDensity+hDensity2+medianIncome, ~hDensity+hDensity2+eduHS, ~age, ~marred,
  ~medianIncome, ~eduHS
  )

names(formulalist) <- as.character(formulalist)


# Function that fits transect models with the formulas in the list above

fit.trans.models <- function(covlist){
  modelist <- vector('list',length=length(covlist))
  nameslist <- names(formulalist)
  for(i in 1:length(covlist)){
    modelist[[i]] <- gdistsamp(
      lambdaformula = covlist[[i]],
      phiformula = ~time,
      pformula = ~temp,
      data = gUmfWithCovs,
      keyfun = 'halfnorm',
      mixture = 'NB'
      )
  }
  names(modelist) <- nameslist
  return(modelist)
}


# Function that fits camera models with the formulas in the list above

fit.cam.models <- function(covlist){
  modelist <- vector('list',length = length(covlist))
  nameslist <- as.character(covlist)
  for(i in 1:length(covlist)){
    modelist[[i]] <- pcount(
      formula = as.formula(
        paste('~dewLow', nameslist[[i]])),
      data = camUmfWithCovs,
      mixture = 'NB',
      K = 50
      )
  }
  names(modelist) <- nameslist
  return(modelist)
}


# Fit the transect models and create an AICc table

transmodels <- fit.trans.models(covlist=formulalist)

transAICTab <- aictab(transmodels, modnames=names(transmodels))


# Fit the camera models and create an AICc table

cammodels <- fit.cam.models(covlist=formulalist)

camAICTab <- aictab(cammodels, modnames=names(cammodels))



# ---------------------------------------------------------------------------------*
# ---- Format transect model results ----
# ---------------------------------------------------------------------------------*

# Find predicted betas for transect models

tAbunds <- vector('list',length=length(transmodels))

for(i in 1:length(tAbunds)){
  tAbunds[[i]] <- predict(transmodels[[i]], type='lambda', appenddata=TRUE)
}


# Weight model betas  and SEs by their AIC weights (transect models)

unsortTransTab <- transAICTab[match(names(transmodels), transAICTab$Modnames),]
tWeight <- unsortTransTab$AICcWt

tV1 <- vector('list',length=length(tAbunds))
tV2 <- vector('list',length=length(tAbunds))


for(i in 1:length(tAbunds)){
  tV1[[i]] <- tAbunds[[i]]$Predicted * tWeight[i]
  tV2[[i]] <- tAbunds[[i]]$SE * tWeight[i]
}

tAbundsWt <- data.frame(matrix(nrow=53,ncol=34))
tSEWt <- data.frame(matrix(nrow=53,ncol=34))
for(i in 1:length(tV1)){
  tAbundsWt[,i] <- tV1[[i]]
}

for(i in 1:length(tV2)){
  tSEWt[,i] <- tV2[[i]]
}

tAbundsWt <- tAbundsWt %>%
  transmute(abund = rowSums(tAbundsWt))

tSEWt <- tSEWt %>%
  transmute(SE = rowSums(tSEWt))

tSites <- catTransect %>%
  select(site) %>%
  distinct

tAbundsWt <- bind_cols(tSites,tAbundsWt,tSEWt)
tAbundsWt[,2:3] <- tAbundsWt[,2:3]/2
colnames(tAbundsWt) <- c('site','abunds','SE')

tAbundsWt <- left_join(tAbundsWt,covs,by='site')


# Examine model-averaged abundance, SE estimates (transect)

tAbundsWt


# Write to a file

write.csv(tAbundsWt, 'data/abundanceTrans.csv', row.names = FALSE)


# ---------------------------------------------------------------------------------*
# ---- Format camera model results ----
# ---------------------------------------------------------------------------------*

# Find model betas for camera models

cAbunds <- vector('list',length=length(cammodels))

for(i in 1:length(cammodels)){
  cAbunds[[i]] <- predict(cammodels[[i]],type='state')
}


# Weight model betas and SEs by their AIC weights, camera models

unsortCamTab <- camAICTab[match(names(cammodels), camAICTab$Modnames),]
cWeight <- unsortCamTab$AICcWt

cV1 <- vector('list',length=length(cAbunds))
cV2 <- vector('list',length=length(cAbunds))


for(i in 1:length(cAbunds)){
  cV1[[i]] <- cAbunds[[i]]$Predicted * cWeight[i]
  cV2[[i]] <- cAbunds[[i]]$SE * cWeight[i]
}

cAbundsWt <- data.frame(matrix(nrow=48,ncol=34))
cSEWt <- data.frame(matrix(nrow=48,ncol=34))
for(i in 1:length(cV1)){
  cAbundsWt[,i] <- cV1[[i]]
}

for(i in 1:length(cV2)){
  cSEWt[,i] <- cV2[[i]]
}

cAbundsWt <- cAbundsWt %>%
  transmute(abund = rowSums(cAbundsWt))

cSEWt <- cSEWt %>%
  transmute(SE = rowSums(cSEWt))

cSites <- catTransect %>%
  select(site) %>%
  filter(!site %in% removeSites) %>%
  distinct

cAbundsWt <- bind_cols(cSites,cAbundsWt,cSEWt)
cAbundsWt[,2:3] <- cAbundsWt[,2:3]/2
colnames(cAbundsWt) <- c('site','abunds','SE')

cAbundsWt <- left_join(cAbundsWt,covs,by='site')


# Examine model-averaged abundance, SE estimates (camera)

cAbundsWt


# Write to a file

write.csv(cAbundsWt, 'data/abundanceCam.csv', row.names = FALSE)


