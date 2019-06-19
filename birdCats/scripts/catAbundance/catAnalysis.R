source('scripts/catAbundance/setup.R')


# prepare data for transect analysis --------------------------------------

# Create an unmarked frame of distance data for the transect counts:

catTransectUmf <- 
  formatDistData(
    data.frame(catTransect), 
    distCol="distance",
    transectNameCol="site", 
    dist.breaks=seq(0, 50, by=5),
    occasionCol='visit'
  )

# Create a dataframe of potential cat abundance covariates by site:

transCovs <- 
  covs %>%
  select(site) %>%
  bind_cols(
    covs %>%
      select(can:eduHS) %>%
      mutate(
        can2 = can^2,
        imp2 = imp^2,
        age2 = age^2,
        medianIncome2 = medianIncome^2,
        hDensity2 = hDensity^2,
        eduHS2 = eduHS^2
      ) %>%
      mutate_all(scaleVar)
  )

# Create a dataframe of potential cat detection covariates by site:

transDetCovs <- 
  catTransect %>%
  select(-c(species:date)) %>%
  distinct %>%
  mutate(time2 = time^2)

# Create a list of wide-form detection covariates:

ySiteCovs <-
  list(
    time = transposeCovariate(transDetCovs, 'time'),
    temp = transposeCovariate(transDetCovs, 'temp'),
    dew = transposeCovariate(transDetCovs, 'dew'),
    time2 = transposeCovariate(transDetCovs, 'time2'),
    doy = transposeCovariate(transDetCovs, 'doy')
  )


# Create unmarkedFrameGDS object for gdistsamp

gUmfWithCovs <-
  unmarkedFrameGDS(
    y = as.matrix(catTransectUmf),
    siteCovs = data.frame(transCovs),
    numPrimary = 6,
    yearlySiteCovs = ySiteCovs,
    survey = 'line',
    dist.breaks =  seq(0, 50, by = 5),
    tlength = rep(200, nrow(catTransectUmf)),
    unitsIn = 'm'
  )

# prepare data for camera analysis --------------------------------------

# Single out transect-only sites:

removeSites <-
  c('OLONMARDC1',
    'WOLFKARDC1',
    'WOLFAMYDC1',
    'GERYERIMD1',
    'MISSEDDC1')

# Remove sites without cameras from covariate data

camCovs <- 
  transCovs %>%
  filter(!site %in% removeSites)


# Create an unmarkedFramePCount object for pcount

camUmfWithCovs <- 
  unmarkedFramePCount(
    umfCam[, -1],
    siteCovs = camCovs,
    obsCovs = camDetCovs
  )

# transect null model fitting ---------------------------------------------

# Rather than building models with every possible combination of phi and p formulas,
# we compare null models, find the best one, use that model's phi and p formulas in 
# the abundance models.

# Function that builds null abundance models with user-selected detection covariates

fit.trans.null.models <-
  function(phi, p) {
    gdistsamp(
      lambdaformula = ~ 1,
      phiformula = phi,
      pformula = p,
      data = gUmfWithCovs,
      keyfun = 'halfnorm',
      mixture = 'NB',
      K = 50
    )
  }

# Formulas for availability:

phiFormulas <-
  c('~time', '~time + time2', '~1')

# Formulas for detection:

pFormulas <-
  c(
    '~doy + temp',
    '~doy + time',
    '~doy + dew',
    '~doy + time + time2',
    '~dew + time',
    '~dew + temp',
    '~dew + time + time2',
    '~time + temp',
    '~time + time2 + temp',
    '~doy',
    '~time',
    '~time + time2',
    '~dew',
    '~temp',
    '~1'
  )


# Create empty list for results:

outList <- 
  vector(
    'list',
    length = length(phiFormulas))


# Fit the models

for (i in 1:length(phiFormulas)) {
  phiList <-
    vector('list', length = length(pFormulas))
  names(phiList) <-
    paste(phiFormulas[i], pFormulas)
  for (j in 1:length(pFormulas)) {
    phiList[[j]] <-
      fit.trans.null.models(phiFormulas[i], pFormulas[j])
  }
  outList[[i]] <- phiList
}

tModels <- 
  unlist(
    outList, 
    recursive = FALSE)


# Make an AIC table to compare detection models

aictab(
  cand.set = tModels,
  modnames = names(tModels))

# camera null model fitting ---------------------------------------------

# Same as transect models, we first compare nulls, then use the best model's 
# p formulas in future abundance models

# Function to fit null abundance models given user-supplied detection formulas:

fit.cam.null.models <- 
  function(p) {
    pcount(
      formula = as.formula(paste(p, '~1')),
      data = camUmfWithCovs,
      mixture = 'NB',
      K = 50
    )
  }


# Detection formulas:

pFormulasCam <-
  c(
    '~doy',
    '~tempHigh',
    '~tempLow',
    '~dewHigh',
    '~doy + tempHigh',
    '~doy + tempLow',
    '~doy + dewHigh',
    '~tempHigh + dewHigh',
    '~tempLow + dewHigh',
    '~1'
  )

# Create an empty vector for results:

cModels <- 
  vector(
    mode = 'list',
    length = length(pFormulasCam)) %>%
  setNames(pFormulasCam)


# Fit the models:

for (i in 1:length(pFormulasCam)) {
  cModels[[i]] <- 
    fit.cam.null.models(pFormulasCam[i])
}


# Make an AIC table to compare camera detection models:

aictab(
  cand.set = cModels,
  modnames = names(cModels))


# abundance model fitting -------------------------------------------------

# Function that fits transect models with user-supplied covariates:

fit.trans.abundance.models <-
  function(formula) {
    gdistsamp(
      lambdaformula = formula,
      phiformula = '~ time',
      pformula = ~ '~doy + dew',
      data = gUmfWithCovs,
      output = 'density',
      unitsOut = 'ha',
      keyfun = 'halfnorm',
      mixture = 'NB'
    )
  }

# Function that fits camera models with user-supplied covariates:

fit.cam.abundance.models <-
  function(formula) {
    pcount(
      formula = as.formula(paste('~dewHigh', as.character(formula))),
      data = camUmfWithCovs,
      mixture = 'NB',
      K = 50
    )
  }


# impermeable surface models ----------------------------------------------

# List of model formulas:

formulaList <-
  c('~ imp',
    '~ imp + imp2',
    '~ 1')

# Create an empty list to store the transect models:

tImpModels <- 
  vector(
    'list',
    length = length(formulaList)) %>%
  setNames(formulaList)

# Fit the transect impermeable models

for (i in 1:length(formulaList)) {
  tImpModels[[i]] <- 
    fit.trans.abundance.models(formulaList[i])
}

# Create an AIC table

tImpAICTab <- 
  aictab(cand.set = tImpModels, modnames=names(tImpModels))


# Create an empty list to store the camera models

cImpModels <- 
  vector(
    'list', 
    length = length(formulaList)) %>%
  setNames(formulaList)


# Fit the camera impermeable models

for (i in 1:length(formulaList)) {
  cImpModels[[i]] <- 
    fit.cam.abundance.models(formulaList[i])
}


# Create an AIC table

cImpAICTab <- 
  aictab(cand.set = cImpModels, modnames=names(cImpModels))



# human demographics models ------------------------------------------------

# List of model formulas

formulaList <-
  c(
    '~1',
    '~hDensity',
    '~hDensity+hDensity2',
    '~medianIncome',
    '~medianIncome+medianIncome2',
    '~eduHS',
    '~eduHS+eduHS2',
    '~medianIncome+eduHS',
    '~hDensity+medianIncome',
    '~hDensity+eduHS',
    '~hDensity+medianIncome+eduHS',
    '~hDensity+hDensity2+medianIncome',
    '~hDensity+hDensity2+eduHS',
    '~hDensity+hDensity2+medianIncome+eduHS',
    '~medianIncome+medianIncome2+hDensity',
    '~medianIncome+medianIncome2+eduHS',
    '~medianIncome+medianIncome2+hDensity+eduHS',
    '~eduHS+eduHS2+hDensity',
    '~eduHS+eduHS2+medianIncome',
    '~eduHS+eduHS2+hDensity+medianIncome',
    '~hDensity+hDensity2+medianIncome+medianIncome2',
    '~hDensity+hDensity2+medianIncome+medianIncome2+eduHS',
    '~hDensity+hDensity2+eduHS+eduHS2',
    '~hDensity+hDensity2+eduHS+eduHS2+medianIncome',
    '~medianIncome+medianIncome2+eduHS+eduHS2',
    '~medianIncome+medianIncome2+eduHS+eduHS2+hDensity',
    '~medianIncome+medianIncome2+eduHS+eduHS2+hDensity+hDensity2'
  )

# Create an empty list to store the transect models:

tDemModels <- 
  vector(
    'list',
    length = length(formulaList)) %>%
  setNames(formulaList)

# Fit the transect impermeable models

for (i in 1:length(formulaList)) {
  tDemModels[[i]] <- 
    fit.trans.abundance.models(formulaList[i])
}

# Create an AIC table

tDemAICTab <- aictab(cand.set = tDemModels, modnames=names(tDemModels))


# Create an empty list to store the camera models

cDemModels <- 
  vector(
    'list', 
    length = length(formulaList)) %>%
  setNames(formulaList)


# Fit the camera impermeable models

for (i in 1:length(formulaList)) {
  cDemModels[[i]] <- 
    fit.cam.abundance.models(formulaList[i])
}


# Create an AIC table

cDemAICTab <- 
  aictab(cand.set = cDemModels, modnames=names(cDemModels))



# Look at the AIC tables

tImpAICTab
tDemAICTab
cImpAICTab
cDemAICTab



# examine model-averaged betasand SEs -------------------------------------------

modavgShrink(tImpModels, 'imp', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(tImpModels, 'imp2', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(tImpModels, 'imp2', parm.type = 'lambda')$Uncond.SE

modavgShrink(cImpModels, 'imp', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(cImpModels, 'imp2', parm.type = 'lambda')$Mod.avg.beta

modavgShrink(tDemModels, 'hDensity', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(tDemModels, 'hDensity2', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(tDemModels, 'eduHS', parm.type = 'lambda')$Mod.avg.beta

modavgShrink(cDemModels, 'hDensity', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(cDemModels, 'hDensity2', parm.type = 'lambda')$Mod.avg.beta
modavgShrink(cDemModels, 'eduHS', parm.type = 'lambda')$Mod.avg.beta



# weight model prediction at each site ----------------------------------------

# Function to weight predicted abundance at each site by model AICc weights

weight.model.predictions <- 
  function(modResults, AICTab, method) {
  
    # Return error for wrong 'method' input
    if (!method %in% c('transect', 'camera')) {
      return("method must be either 'transect' or 'camera'")
    }
    
    # Get appropriate site list depending on detection method
    if (method == 'transect'){
      sites <- covs %>% select(site)
    } else {
      sites <- camCovs %>% select(site)
    }
    
    # Resort AIC table results to match order of modResults
    unsortAICWt <- AICTab[match(names(modResults), AICTab$Modnames),]$AICcWt
    
    # Make empty list to store weighted model predictions
    weighted <- vector('list', length = length(modResults))
    names(weighted) <- names(modResults)
    
    # Multiply model predictions by their AICc weights
    for (i in 1:length(modResults)) {
      weighted[[i]] <- modResults[[i]] * unsortAICWt[i]
    }
    
    # Add up weighted predictions
    for (i in 1:ncol(weighted[[1]])) {
      for (j in 2:length(weighted)) {
        weighted[[1]][i] <- weighted[[1]][i] + weighted[[j]][i]
      }
    }
    
    # Set up output dataframe
    output <- bind_cols(sites, weighted[[1]])
    colnames(output) <- c('site', 'abund', 'SE', 'lower', 'upper')
    
    # Add covariates to output dataframe
    output <- left_join(output, covs, by = 'site')
    
    return(output)
}


# Find predicted abundances for transect models

tImpUnweightAbund <- 
  lapply(
    tImpModels, 
    predict, 
    type = 'lambda', 
    appenddata = TRUE)

tDemUnweightAbund <- 
  lapply(
    tDemModels, 
    predict, 
    type = 'lambda', 
    appenddata = TRUE)

cImpUnweightAbund <- 
  lapply(
    cImpModels, 
    predict, 
    type = 'state', 
    appenddata = TRUE)

cDemUnweightAbund <- 
  lapply(
    cDemModels, 
    predict, 
    type = 'state', 
    appenddata = TRUE)



# Weight predicted abundance by AIcC weights

tImpWeightAbund <- 
  weight.model.predictions(
    tImpUnweightAbund, 
    tImpAICtab, 
    'transect')

tDemWeightAbund <- 
  weight.model.predictions(
    tDemUnweightAbund, 
    tDemAICtab, 
    'transect')

cImpWeightAbund <- 
  weight.model.predictions(
    cImpUnweightAbund, 
    cImpAICtab, 
    'camera')

cDemWeightAbund <- 
  weight.model.predictions(
    cDemUnweightAbund, 
    cDemAICtab, 
    'camera')


