# =================================================================================*
# ---- SET-UP ----
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

smartLibrary(c('unmarked', 'dplyr', 'tidyr', 'camtrapR', 'ggplot2', 'AICcmodavg','MuMIn'))


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
    eduC = scale(eduC)[,1],
    marred = scale(marred)[,1])

gUmfWithCovs@yearlySiteCovs <- gUmfWithCovs@yearlySiteCovs %>%
  mutate(
    time2 = time^2,
    time = scale(time)[,1],
    time2 = scale(time2)[,1],
    temp = scale(temp)[,1],
    dew = scale(dew)[,1])


# ---------------------------------------------------------------------------------*
# ---- Transect model fitting ----
# ---------------------------------------------------------------------------------*


# Null model

null.model <-  gdistsamp(
  lambdaformula = ~1,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# Global model

global.model <- gdistsamp(
  lambdaformula = ~imp+imp2+age+age2+eduC+marred+medianIncome+medianIncome2,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")



# Single covariates

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


# # median income
# 
# income.model <- gdistsamp(
#   lambdaformula = ~ medianIncome,
#   phiformula = ~ time,
#   pformula = ~ dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")
# 
# 
# # age
# 
# age.model <- gdistsamp(
#   lambdaformula = ~age,
#   phiformula = ~time,
#   pformula = ~dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")


# human density

density.model <- gdistsamp(
  lambdaformula = ~hDensity,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# # married
# 
# marred.model <- gdistsamp(
#   lambdaformula = ~marred,
#   phiformula = ~time,
#   pformula = ~dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")
# 
# 
# # HS education
# 
# eduHS.model <- gdistsamp(
#   lambdaformula = ~eduHS,
#   phiformula = ~time,
#   pformula = ~dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")
# 
# 
# # C education
# 
# eduC.model <- gdistsamp(
#   lambdaformula = ~eduC,
#   phiformula = ~time,
#   pformula = ~dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")



# Additive

# # canopy + medianIncome
# 
# can.income.model <- gdistsamp(
#   lambdaformula = ~ can+medianIncome,
#   phiformula = ~ time,
#   pformula = ~ dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")
# 
# 
# # imp + medianIncome
# 
# imp.income.model <- gdistsamp(
#   lambdaformula = ~ imp+medianIncome,
#   phiformula = ~ time,
#   pformula = ~ dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm",  mixture="NB")
# 

# Quadratic

# # income^2
# 
# income2.model <- gdistsamp(
#   lambdaformula = ~ medianIncome+medianIncome2,
#   phiformula = ~ time,
#   pformula = ~ dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")
# 


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


# impervious^2 + median income

imp2.income.model <- gdistsamp(
  lambdaformula = ~imp+imp2+medianIncome,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# # canopy^2 + medianIncome
# 
# can2.income.model <- gdistsamp(
#   lambdaformula = ~ can+can2+medianIncome,
#   phiformula = ~ time,
#   pformula = ~ dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")


# hDensity^2

density2.model <- gdistsamp(
  lambdaformula = ~hDensity+hDensity2,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# impervious^2 + eduHS

imp2.eduHS.model <- gdistsamp(
  lambdaformula = ~imp+imp2+eduHS,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# impervious^2 + eduC <- this is the best model

imp2.eduC.model <- gdistsamp(
  lambdaformula = ~imp+imp2+eduC,
  phiformula = ~time,
  pformula = ~dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# # impervious^2 + income + C education
# 
# imp2.income.eduC.model <- gdistsamp(
#   lambdaformula = ~imp+imp2+eduC+medianIncome,
#   phiformula = ~time,
#   pformula = ~dew*temp+time,
#   data = gUmfWithCovs,
#   keyfun = "halfnorm", mixture="NB")

# impervious^2 + income^2

imp2.income2.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+medianIncome+medianIncome2,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# impervious^2 + age

imp2.age.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+age,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")

# impervious^2 + age^2

imp2.age2.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+age+age2,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")


# impervious^2 + marred

imp2.marred.model <- gdistsamp(
  lambdaformula = ~ imp+imp2+marred,
  phiformula = ~ time,
  pformula = ~ dew*temp+time,
  data = gUmfWithCovs,
  keyfun = "halfnorm", mixture="NB")



# -------------------*
# ---- Abundance ----
# -------------------*


# Model selection tables

modelList1 <- list(imp.model,can.model,density.model,
                   imp2.model,can2.model,density2.model)
names1 <- c('imp','can','hDensity','imp2','can2','hDensity2')

table1 <- aictab(cand.set=modelList1,modnames=names1)


modelList2 <- list(imp2.model,imp2.age.model,imp2.income.model,imp2.income2.model,
                   imp2.marred.model,imp2.age2.model,imp2.eduHS.model,imp2.eduC.model)
names2 <- c('imp2','imp2.age','imp2.income','imp2.income2','imp2.marred','imp2.age2',
            'imp2.eduHS','imp2.eduC')

table2 <- aictab(cand.set=modelList2,modnames=names2)

tAvg <- model.avg(modelList2)
tWeight <- tAvg$msTable$weight
seImp <- tAvg$coefArray[,2,2]
seImp2 <- tAvg$coefArray[,2,3]

seImpWt <- seImp*tWeight
seImp2Wt <- seImp2*tWeight

sum(seImpWt)
sum(seImp2Wt)



# Predicted abund values from best model

transAbund <- predict(imp2.eduHS.model, type="lambda", appenddata = TRUE)




# ----------------------------------------------------------------*
# ------ Visualize variables ------
# ----------------------------------------------------------------*



plot(covs$imp, transAbund$Predicted)
plot(covs$eduHS, transAbund$Predicted)




# =============================================================================*
# ---- CAMERA DATA ----
# =============================================================================*

# Create detection history for each site:

umfCam <- read.csv('data/catCamDetection.csv') %>%
  data.frame


# Get abundance covariates for the camera-only sites
# Square and scale variables

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1',
                 'GERYERIMD1', 'MISSEDDC1')

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
    eduC = scale(eduC)[,1],
    marred = scale(marred)[,1])


# Get scaled detection covariates

camDetCovs <- read.csv('data/camDetCovs.csv') %>%
  filter(!site %in% removeSites) %>%
  select(site, day, tempHigh, tempLow, dewLow) %>%
  mutate(
    tempHigh = scale(tempHigh)[,1],
    tempLow = scale(tempLow)[,1],
    dewLow = scale(dewLow)[,1])


# Create an unmarkedFramePCount object for pcount

camUmfWithCovs <- unmarkedFramePCount(
  umfCam[,-1],
  siteCovs = camCovs,
  obsCovs = camDetCovs)


# ------------------------------------------------------------------------------*
# ---- Cam model fitting ----
#-------------------------------------------------------------------------------*


#Null

null.cmodel <- pcount(
  formula = ~dewLow ~1, 
  data= camUmfWithCovs,
  mixture = 'NB', K = 50)


# Global

global.cmodel <- pcount(
  formula = ~dewLow
    ~can + can2 + medianIncome + medianIncome2 + age + age2 + marred + eduC + eduHS,
  data= camUmfWithCovs, mixture = 'NB', K = 50)


# Single abundance covariates

imp.cmodel <- pcount(
  formula = ~dewLow ~imp, 
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can.cmodel <- pcount(
  formula = ~dewLow ~can, 
  data = camUmfWithCovs, mixture = 'NB', K = 50)

density.cmodel <- pcount(
  formula = ~dewLow ~hDensity, 
  data = camUmfWithCovs, mixture = 'NB', K = 50)

# income.cmodel <- pcount(
#   formula = ~dewLow ~medianIncome,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# eduC.cmodel <- pcount(
#   formula = ~dewLow ~eduC,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# eduHS.cmodel <- pcount(
#   formula = ~dewLow ~eduHS,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# age.cmodel <- pcount(
#   formula = ~dewLow ~age,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# marred.cmodel <- pcount(
#   formula = ~dewLow ~marred,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)



# Additive

# can.income.cmodel <- pcount(
#   formula = ~dewLow ~can + medianIncome,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# imp.income.cmodel <- pcount(
#   formula = ~dewLow ~imp + medianIncome,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)



# Quadratic

imp2.cmodel <- pcount(
  formula = ~dewLow ~imp + imp2,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.cmodel <- pcount(
  formula = ~dewLow ~can + can2,
  data = camUmfWithCovs, mixture= 'NB', K = 50)

density2.cmodel <- pcount(
  formula = ~dewLow ~hDensity + hDensity2,
  data = camUmfWithCovs, mixture= 'NB', K = 50)

# age2.cmodel <- pcount(
#   formula = ~dewLow ~age+age2,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)
# 
# imp2.income.cmodel <- pcount(
#   formula = ~dewLow ~imp+imp2+medianIncome,
#   data = camUmfWithCovs, mixture='NB', K = 50)

can2.income.cmodel <- pcount(
  formula = ~dewLow ~can + can2 + medianIncome,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

# imp2.eduC.cmodel <- pcount(
#   formula = ~dewLow ~imp + imp2 + eduC,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.eduC.cmodel <- pcount(
  formula = ~dewLow ~can + can2 + eduC,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

# imp2.eduHS.cmodel <- pcount(
#   formula = ~dewLow ~imp + imp2 + eduHS,
#   data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.eduHS.cmodel <- pcount(
  formula = ~dewLow ~can + can2 + eduHS,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.marred.cmodel <- pcount(
  formula = ~dewLow ~can+can2+marred,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.age.cmodel <- pcount(
  formula = ~dewLow ~can+can2+age,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.age2.cmodel <- pcount(
  formula = ~dewLow ~can+can2+age+age2,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.income2.cmodel <- pcount(
  formula = ~dewLow ~can+can2+medianIncome+medianIncome2,
  data = camUmfWithCovs, mixture = 'NB', K = 50)

can2.income.eduC.cmodel <- pcount(
  formula = ~dewLow ~can+can2+medianIncome+eduC,
  data = camUmfWithCovs, mixture = 'NB', K = 50)



# -------------------*
# ---- Abundance ----
# -------------------*

# Model selection tables

cmodelList1 <- list(imp.cmodel,can.cmodel,density.cmodel,
                   imp2.cmodel,can2.cmodel,density2.cmodel)
cnames1 <- c('imp','can','hDensity','imp2','can2','hDensity2')

ctable1 <- aictab(cand.set=cmodelList1,modnames=cnames1)


cmodelList2 <- list(can2.cmodel,can2.age.cmodel,can2.income.cmodel,
                    can2.income2.cmodel,can2.marred.cmodel,can2.age2.cmodel,
                    can2.eduHS.cmodel,can2.eduC.cmodel)
cnames2 <- c('can2','can2.age','can2.income','can2.income2','can2.marred',
             'can2.age2','can2.eduHS','can2.eduC')

ctable2 <- aictab(cand.set=cmodelList2,modnames=cnames2)

cAvg <- model.avg(cmodelList2)
cWeight <- cAvg$msTable$weight
seCan <- cAvg$coefArray[,2,2]
seCan2 <- cAvg$coefArray[,2,3]

seCanWt <- seCan*tWeight
seCan2Wt <- seCan2*tWeight

sum(seCanWt)
sum(seCan2Wt)


# Predicted abundance based on best model

camAbund <- predict(can2.eduHS.cmodel, type = 'state')


# -------------------------------------------------*
# ------ Visualize variables ------
# -------------------------------------------------*

camSites <- covs %>% filter(!site %in% removeSites)

plot(camSites$eduHS, camAbund$Predicted)
plot(camSites$can, camAbund$Predicted)


# Everything below this point is commented out


# # ---------------------------------------*
# # ----- Combine abundance estimates -----
# # ---------------------------------------*
# 
# # Create dataframe of camera sites and abundances
# 
# camSiteAbund <- sitesWithCovs %>%
#   filter(!site %in% removeSites) %>%
#   mutate(cCats = camAbund$Predicted) %>%
#   select(site, cCats)
# 
# 
# # Read in minimum individual data
# 
# minInd <- read.csv('data/catMinIndividuals.csv') %>%
#   select(site, mCats)
# 
# 
# # Create dataframe of transect sites and predicted abundances
# 
# siteAbund <- sitesWithCovs %>%
#   select(site) %>%
#   mutate(tCats = transAbund$Predicted) %>%
#   left_join(camSiteAbund, by = 'site') %>%
#   left_join(minInd, by = 'site')
# 
# 
# # Write to csv file for use in survival analysis
# 
# write.csv(siteAbund, 'data/catAbund.csv', row.names = FALSE)
# 
# 




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


# # -------------*
# # ---- PLOT ----
# # -------------*
# 
# sitesImp <- select(sitesWithCovs, site, imp)
# 
# indsImpData <- indsPerSite %>%
#   left_join(sitesImp, by = 'site')
# 
# indsImpData$site <- as.factor(indsImpData$site)
# 
# ggplot(data = indsImpData, aes(x = imp, y = total))+
#   geom_point(size = 2)+
#   theme(panel.grid = element_blank(), 
#         axis.line.x = element_line(linetype='solid', color='black'),
#         axis.line.y = element_line(linetype='solid', color='black'))
#   
#   
#   
  
