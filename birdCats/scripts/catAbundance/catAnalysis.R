# =================================================================================*
# --------------------------------- Set-up ----------------------------------------
# =================================================================================*

# ------------------------*
# ------ Functions -------
# ------------------------*

# Function searches packages in installed package list, installs them if they are 
# not present, and loads the library:

smartLibrary <- function(packageVector) {
  for (i in 1:length(packageVector)) {
    package <- packageVector[i]
    if (!package %in% rownames(installed.packages())) {
      install.packages(packageVector[i],
                       repos = "http://cran.rstudio.com/",
                       dependencies = TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}


# Function that scales variables

scaleVar <- function(var) {
  (var - mean(var, na.rm = TRUE))/sd(var, na.rm = TRUE)
}


# Function that scales and transposes detection covariates to wide format

transposeCovariate <- function(data, covariate){
  transMat <- data %>%
    select(site, visit, cov = covariate) %>%
    mutate(cov = scaleVar(cov)) %>%
    spread(visit, cov) %>%
    select(-site) %>%
    as.matrix()
  colnames(transMat) <- NULL
  return(transMat)
}

# ------------------------*
# --- Load libraries -----
# ------------------------*

smartLibrary(
  c('unmarked', 'dplyr', 'tidyr', 'camtrapR', 'ggplot2', 'AICcmodavg','MuMIn')
)

# ------------------------*
# ---- Load the data -----
# ------------------------*

options(stringsAsFactors = F)


# Site visit data: time and weather

catSiteActivity <- read.csv('data/catDataActivity.csv') %>%
  tbl_df


# Camera sampling data

catCam <- read.csv('data/catDataCamera.csv') %>%
  tbl_df %>%
  filter(!is.na(species), !is.na(cameraID))


# Total Neighborhood Nestwatch site list

catSites <- read.csv('data/catSiteData.csv') %>%
  tbl_df %>%
  select(site)


# Transect sampling data

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


# Education data by site

eduData <- read.csv('data/eduHS.csv')


# Covariate data by site

covs <- read.csv('data/covariateData.csv') %>%
  select(-c(eduHS, eduC)) %>%
  left_join(eduData) %>%
  tbl_df %>%
  arrange(site)


# Camera detection history of cats for each site:

umfCam <- read.csv('data/catCamDetection.csv')


# Scaled camera detection covariates

camDetCovs <- read.csv('data/camDetCovs.csv') %>%
  mutate(
    tempHigh = scaleVar(tempHigh),
    tempLow = scaleVar(tempLow),
    dewLow = scaleVar(dewLow)
  )


# =================================================================================*
# --------------------------- Set up transect analysis ----------------------------
# =================================================================================*

# Create an unmarked frame of distance data for the transect counts

catTransectUmf <- formatDistData(
  data.frame(catTransect), 
  distCol="distance",
  transectNameCol="site", 
  dist.breaks=seq(0, 50, by=5),
  occasionCol='visit'
)


# Create a dataframe of potential cat abundance covariates by site

transCovs <- covs %>%
  select(site) %>%
  bind_cols(
    covs %>%
      select(can:marred) %>%
      mutate(
        can2 = can^2,
        imp2 = imp^2,
        age2 = age^2,
        medianIncome = medianIncome^2,
        hDensity2 = hDensity^2
      ) %>%
      mutate_all(scaleVar)
  )


# Create a dataframe of potential cat detection covariates by site

transDetCovs <- catTransect %>%
  select(-c(species:date)) %>%
  distinct %>%
  mutate(time2 = time^2)


# Create a list of wide-form detection covariates

ySiteCovs <- list(
  time = transposeCovariate(transDetCovs, 'time'),
  temp = transposeCovariate(transDetCovs, 'temp'),
  dew = transposeCovariate(transDetCovs, 'dew'),
  time2 = transposeCovariate(transDetCovs, 'time2')
)


# Create unmarkedFrameGDS object for gdistsamp

gUmfWithCovs <- unmarkedFrameGDS(
  y =              as.matrix(catTransectUmf),
  siteCovs =       data.frame(transCovs),
  numPrimary =     6,
  yearlySiteCovs = ySiteCovs,
  survey =         'line',
  dist.breaks=     seq(0, 50, by=5),
  tlength =        rep(200, nrow(catTransectUmf)),
  unitsIn =        'm'
)


# =================================================================================*
# ------------------------- Transect null model fitting ---------------------------
# =================================================================================*

# Rather than building models with every possible combination of phi and p
# formulas, we compare null models, find the best one, use that model's phi and
# p formulas in the abundance models.


# Function that builds null abundance models with user-selected
# detection covariates

nullModelFit <- function(phi, p) {
  gdistsamp(
    lambdaformula = ~1,
    phiformula = phi,
    pformula = p,
    data = gUmfWithCovs,
    keyfun = 'halfnorm',
    mixture = 'NB'
  )
}


# Formulas for availability

phiFormulas <- c('~time', '~1')


# Formulas for detection

pFormulas <- c(
  '~dew + temp + time',
#  '~dew + time', (NaNs produced)
  '~temp + time',
  '~dew + temp',
  '~time',
#  '~dew', (NaNs produced)
  '~temp',
  '~1'
)


# Create empty list for results

outList <- vector('list', length = length(phiFormulas))


# Fit the models

for (i in 1:length(phiFormulas)) {
  phiList <- vector('list', length = length(pFormulas))
  names(phiList) <- paste(phiFormulas[i], pFormulas)
  for (j in 1:length(pFormulas)) {
    phiList[[j]] <- nullModelFit(phiFormulas[i], pFormulas[j])
  }
  outList[[i]] <- phiList
}

modelList <- unlist(outList, recursive = FALSE)


# Make an AIC table to compare detection models

aictab(cand.set = modelList, modnames = names(modelList))



# =================================================================================*
# ---------------------------- Set up camera analysis -----------------------------
# =================================================================================*

# Single out transect-only sites

removeSites <- c('OLONMARDC1','WOLFKARDC1', 'WOLFAMYDC1', 'GERYERIMD1', 'MISSEDDC1')


# Remove sites without cameras from covariate data

camCovs <- transCovs %>%
  filter(!site %in% removeSites)


# Create an unmarkedFramePCount object for pcount

camUmfWithCovs <- unmarkedFramePCount(
  umfCam[, -1],
  siteCovs = camCovs,
  obsCovs = camDetCovs
)


# =================================================================================*
# ------------------------- Camera null model fitting -----------------------------
# =================================================================================*

# Same as transect models, we first compare nulls, then use the best model's p 
# formulas in future abundance models


# Function to fit null abundance models given user-supplied detection formulas

nullCamModelFit <- function(p) {
  pcount(
    formula = as.formula(paste(p, '~1')),
    data = camUmfWithCovs,
    mixture = 'NB',
    K = 50
  )
}


# Detection formulas

pFormulasCam <- c(
  '~tempHigh',
  '~tempLow',
  '~dewLow',
  '~tempHigh + dewLow',
  '~tempLow + dewLow',
  '~tempHigh * dewLow',
  '~tempLow * dewLow',
  '~1'
)


# Create an empty vector for results

modelListCam <- vector(mode = 'list', length = length(pFormulasCam))
names(modelListCam) <- pFormulasCam


# Fit the models

for (i in 1:length(pFormulasCam)) {
  modelListCam[[i]] <- nullCamModelFit(pFormulasCam[i])
}


# Make an AIC table to compare camera detection models

aictab(cand.set = modelListCam, modnames = names(modelListCam))


# =================================================================================*
# ---------------------------- Abundance model fitting ----------------------------
# =================================================================================*

# Function that fits transect models with user-supplied covariates

fit.trans.models <- function(formula) {
  gdistsamp(
    lambdaformula = formula,
    phiformula = ~time,
    pformula = ~1,
    data = gUmfWithCovs,
    keyfun = 'halfnorm',
    mixture = 'NB'
  )	
}


# Function that fits camera models with user-supplied covariates

fit.cam.models <- function(formula) {
  pcount(
    formula = as.formula(paste('~dewLow', as.character(formula))),
    data = camUmfWithCovs,
    mixture = 'NB',
    K = 50
  )
}


# List of model formulas

formulaList <- c(
  ~imp,
  ~imp+imp2, 
  ~imp+marred,
  ~imp+medianIncome,
  ~imp+imp2+marred,
  ~imp+imp2+medianIncome, 
  ~marred,
  ~medianIncome
)


# Create an empty list to store the transect models

tModels <- vector('list', length = length(formulaList))
names(tModels) <- as.character(formulaList)


# Fit the transect models

for (i in 1:length(formulaList)) {
  tModels[[i]] <- fit.trans.models(formulaList[i])
}


# Create an AIC table

tAICTab <- aictab(cand.set = tModels, modnames=names(tModels))


# Create an empty list to store the camera models

cModels <- vector('list', length = length(formulaList))
names(cModels) <- as.character(formulaList)


# Fit the camera models

for (i in 1:length(formulaList)) {
  cModels[[i]] <- fit.cam.models(formulaList[i])
}


# Create an AIC table

cAICTab <- aictab(cand.set = cModels, modnames=names(cModels))



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

tAbundsWt <- data.frame(matrix(nrow=nrow(covs), ncol=length(tAbunds)))
tSEWt <- data.frame(matrix(nrow=nrow(covs), ncol=length(tAbunds)))
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

cAbundsWt <- data.frame(matrix(nrow=nrow(camCovs), ncol=length(cAbunds)))
cSEWt <- data.frame(matrix(nrow=nrow(camCovs), ncol=length(cAbunds)))
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

cSites <- camCovs %>%
  select(site)

cAbundsWt <- bind_cols(cSites,cAbundsWt,cSEWt)
cAbundsWt[,2:3] <- cAbundsWt[,2:3]/2
colnames(cAbundsWt) <- c('site','abunds','SE')

cAbundsWt <- left_join(cAbundsWt,covs,by='site')


# Examine model-averaged abundance, SE estimates (camera)

cAbundsWt


# Write to a file

write.csv(cAbundsWt, 'data/abundanceCam.csv', row.names = FALSE)


