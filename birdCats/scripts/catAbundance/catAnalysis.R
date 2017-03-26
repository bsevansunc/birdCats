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

covs <- read.csv('data/covariateData.csv') %>%
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

# Scale variables

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
    eduC = scale(eduC)[,1],
    marred = scale(marred)[,1])

gUmfWithCovs@yearlySiteCovs <- gUmfWithCovs@yearlySiteCovs %>%
  mutate(
    time = scale(time)[,1],
    temp = scale(temp)[,1],
    dew = scale(dew)[,1])



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
# ---- Transect model fitting ----
# ---------------------------------------------------------------------------------*

# gdistsamp

# Null model

null.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Global model

global.model <- gdistsamp(
  lambdaformula = ~imp+imp2+age+marred+medianIncome+medianIncome2,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# impervious^2 + median income

imp2.income.model <- gdistsamp(
  lambdaformula = ~imp+imp2+medianIncome,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# impervious

imp.model <- gdistsamp(
  lambdaformula = ~imp,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# canopy

can.model <- gdistsamp(
  lambdaformula = ~can,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# median income

income.model <- gdistsamp(
  lambdaformula = ~ medianIncome,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# income2

income2.model <- gdistsamp(
  lambdaformula = ~ medianIncome+medianIncome2,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")



# impervious^2

imp2.model <- gdistsamp(
  lambdaformula = ~ imp + imp2,
  phiformula = ~ time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# canopy^2

can2.model <- gdistsamp(
  lambdaformula = ~can+can2,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs, mixture="NB")


# canopy^2 + medianIncome

can2.income.model <- gdistsamp(
  lambdaformula = ~ can+can2+medianIncome,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# canopy + medianIncome

can.income.model <- gdistsamp(
  lambdaformula = ~ can+medianIncome,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# imp + medianIncome

imp.income.model <- gdistsamp(
  lambdaformula = ~ imp+medianIncome,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm",  mixture="NB")


# age

age.model <- gdistsamp(
  lambdaformula = ~age,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# human density

density.model <- gdistsamp(
  lambdaformula = ~hDensity,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# hDensity2

density2.model <- gdistsamp(
  lambdaformula = ~hDensity+hDensity2,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# married

marred.model <- gdistsamp(
  lambdaformula = ~marred,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# HS education

eduHS.model <- gdistsamp(
  lambdaformula = ~eduHS,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# C education

eduC.model <- gdistsamp(
  lambdaformula = ~eduC,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# imp2 eduHS

imp2.eduHS.model <- gdistsamp(
  lambdaformula = ~imp+imp2+eduHS,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# imp2 eduC

imp2.eduC.model <- gdistsamp(
  lambdaformula = ~imp+imp2+eduC,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# imp2 + income + C education

imp2.income.eduC.model <- gdistsamp(
  lambdaformula = ~imp+imp2+eduC+medianIncome,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# imp2 + income2

imp2.income2.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+medianIncome+medianIncome2,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# imp2 + age

imp2.age.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+age,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# imp2 + age2

imp2.age2.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+age+age2,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# imp2 + marred

imp2.marred.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+marred,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Model selection table

modelList1 <- list(imp.model,can.model,density.model,
                   imp2.model,can2.model,density2.model)
names1 <- c('imp','can','hDensity','imp2','can2','hDensity2')

table1 <- aictab(cand.set=modelList1,modnames=names1)

# Model selection table

modelList2 <- list(imp2.model,imp2.age.model,imp2.income.model,imp2.income2.model,
                   imp2.marred.model,imp2.age2.model,imp2.eduHS.model,imp2.eduC.model)
names2 <- c('imp2','imp2.age','imp2.income','imp2.income2','imp2.marred','imp2.age2',
            'imp2.eduHS','imp2.eduC')

table2 <- aictab(cand.set=modelList2,modnames=names2)

# -------------------*
# ---- Abundance ----
# -------------------*

transAbund <- predict(imp2.income2.model, type="lambda", appenddata = TRUE)



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
  select(site, imp, can, medianIncome, eduC,hDensity)

transSiteAbund <- cbind.data.frame(transSites, transAbund$Predicted)
can <- transSiteAbund$can
imp <- transSiteAbund$imp
inc <- transSiteAbund$medianIncome
abund <- transAbund$Predicted
edu <- transSiteAbund$eduC
dens <- transSiteAbund$hDensity

plot(imp, abund)
plot(edu, abund)
plot(inc, abund)
plot(dens, abund)




# =================================================================================*
# ---- CAMERA DATA ----
# =================================================================================*

# Create detection history for each site:

umfCam <- read.csv('data/catCamDetection.csv') %>%
  data.frame


# Get abundance covariates for the camera-only sites

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1','GERYERIMD1', 'MISSEDDC1')

camCovs <- covs %>%
  filter(!site %in% removeSites) %>%
  data.frame


# Get detection covariates

camDetCovs <- read.csv('data/camDetCovs.csv') %>%
  filter(!site %in% removeSites) %>%
  select(site, day, tempHigh, tempLow, dewLow)


# Create and unmarkedFramePCount object for pcount

camUmfWithCovs <- unmarkedFramePCount(
  umfCam[,-1],
  siteCovs = camCovs,
  obsCovs = camDetCovs
  )


# ---------------------------------------------------------------------------------*
# ---- Cam model fitting ----
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


# -------------------*
# ---- Abundance ----
# -------------------*

camAbund <- predict(camCan2, type = 'state')


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


# ---------------------------------------*
# ----- Combine abundance estimates -----
# ---------------------------------------*

# Create dataframe of camera sites and abundances

camSiteAbund <- sitesWithCovs %>%
  filter(!site %in% removeSites) %>%
  mutate(cCats = camAbund$Predicted) %>%
  select(site, cCats)


# Read in minimum individual data

minInd <- read.csv('data/catMinIndividuals.csv') %>%
  select(site, mCats)


# Create dataframe of transect sites and predicted abundances

siteAbund <- sitesWithCovs %>%
  select(site) %>%
  mutate(tCats = transAbund$Predicted) %>%
  left_join(camSiteAbund, by = 'site') %>%
  left_join(minInd, by = 'site')


# Write to csv file for use in survival analysis

write.csv(siteAbund, 'data/catAbund.csv', row.names = FALSE)






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

# ----------------------------------------------------*
# ---- Minimum number of individuals ----
# ----------------------------------------------------*

na.zero <- function(x){
  x[is.na(x)] <- 0
  return(x)
}

transInds <- catTransect %>%
  filter(species == 'cat' & count == 1) %>%
  select(site, notes) %>%
  unique

transInds$site <- as.character(transInds$site)
transInds$notes <- as.character(transInds$notes)


camInds <- catCam %>%
  filter(species == 'cat' & count == 1) %>%
  select(site, note) %>%
  rename(notes = note) %>%
  unique

camInds$site <- as.character(camInds$site)
camInds$notes <- as.character(camInds$notes)


totalInds <- union(camInds, transInds) %>%
  arrange(by = site) %>%
  filter(notes != '')

sites <- catSiteActivity %>%
  select(site) %>%
  unique

totalIndsPerSite <- data.frame(table(totalInds$site)) %>%
  rename(site = Var1, inds = Freq) %>%
  full_join(sites, by = 'site') %>%
  na.zero() %>%
  arrange(by = site)

camIndsPerSite <- data.frame(table(camInds$site)) %>%
  rename(site = Var1, inds = Freq) %>%
  full_join(sites, by = 'site') %>%
  na.zero() %>%
  arrange(by = site)

transIndsPerSite <- data.frame(table(transInds$site)) %>%
  rename(site = Var1, inds = Freq) %>%
  full_join(sites, by = 'site') %>%
  na.zero() %>%
  arrange(by = site)

indsPerSite <- left_join(totalIndsPerSite, camIndsPerSite, by = 'site') %>%
  left_join(transIndsPerSite, by = 'site') %>%
  rename(total = inds.x, camera = inds.y, transect = inds)

# -------------*
# ---- PLOT ----
# -------------*

sitesImp <- select(sitesWithCovs, site, imp)

indsImpData <- indsPerSite %>%
  left_join(sitesImp, by = 'site')

indsImpData$site <- as.factor(indsImpData$site)

ggplot(data = indsImpData, aes(x = imp, y = total))+
  geom_point(size = 2)+
  theme(panel.grid = element_blank(), 
        axis.line.x = element_line(linetype='solid', color='black'),
        axis.line.y = element_line(linetype='solid', color='black'))
  
  
  
  
