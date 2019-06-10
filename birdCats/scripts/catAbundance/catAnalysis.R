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

# Rather than building models with every possible combination of phi and p formulas, we compare null models, find the best one, use that model's phi and p formulas in the abundance models.

# Function that builds null abundance models with user-selected detection covariates

nullModelFit <-
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
    '~doy + + time + time2',
    '~dew + time',
    '~dew + temp',
    '~dew + time + time2',
    '~time + temp',
    '~time + time2 + temp',
    '~doy',
    '~time',
    '~time + time2',
    '~dew',
    # This term is the one that produces NaNs
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
      nullModelFit(phiFormulas[i], pFormulas[j])
  }
  outList[[i]] <- phiList
}

modelList <- 
  unlist(
    outList, 
    recursive = FALSE)


# Make an AIC table to compare detection models

aictab(
  cand.set = modelList,
  modnames = names(modelList))

# camera null model fitting ---------------------------------------------

# Same as transect models, we first compare nulls, then use the best model's p formulas in future abundance models

# Function to fit null abundance models given user-supplied detection formulas:

nullCamModelFit <- 
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

modelListCam <- 
  vector(
    mode = 'list',
    length = length(pFormulasCam))

names(modelListCam) <- 
  pFormulasCam


# Fit the models:

for (i in 1:length(pFormulasCam)) {
  modelListCam[[i]] <- 
    nullCamModelFit(pFormulasCam[i])
}


# Make an AIC table to compare camera detection models:

aictab(
  cand.set = modelListCam,
  modnames = names(modelListCam))


# abundance model fitting -------------------------------------------------

# Function that fits transect models with user-supplied covariates:

fit.trans.models <-
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

fit.cam.models <-
  function(formula) {
    pcount(
      formula = as.formula(paste('~dewHigh', as.character(formula))),
      data = camUmfWithCovs,
      mixture = 'NB',
      K = 50
    )
  }

# List of model formulas:

formulaList <-
  c('~ imp',
   ' ~ imp + imp2',
    '~ 1')

# Create an empty list to store the transect models:

tModels <- 
  vector(
    'list',
    length = length(formulaList)) %>%
  setNames(formulaList)

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

cAICTab <- 
  aictab(cand.set = cModels, modnames=names(cModels))



# --------------------------------------------------------------------------------*
# --------- Combine models for detection and abundance -----------
# --------------------------------------------------------------------------------*

fit.all.trans.models <- function(lformula, phformula, pformula) {
  gdistsamp(
    lambdaformula = lformula,
    phiformula = phformula,
    pformula = pformula,
    data = gUmfWithCovs,
    output = 'density',
    unitsOut = 'ha',
    keyfun = 'halfnorm',
    mixture = 'NB'
  )	
}

fit.all.cam.models <- function(p.formula, l.formula) {
  pcount(
    formula = as.formula(paste(p.formula, l.formula)),
    data = camUmfWithCovs,
    mixture = 'NB',
    K = 50
  )
}


# Model parameters for both transects and cameras

l.formulas.imp <-
  c('~1',
    '~imp',
    '~imp+imp2')

l.formulas.dem <-
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

ph.formulas <-
  c('~1',
    '~time',
    '~time+time2')

p.formulas.trans <-
  c('~1',
    '~doy',
    '~time',
    '~temp',
    '~dew')

p.formulas.cam <-
  c('~1',
    '~doy',
    '~tempHigh',
    '~tempLow',
    '~dewHigh')


# Fit the models for transects/impervious

outList.imp.trans <- vector('list', length = length(ph.formulas))

for (i in 1:length(ph.formulas)) {
  phiList <- vector('list', length = length(p.formulas.trans))
  for (j in 1:length(p.formulas.trans)) {
    pList <- vector('list', length = length(l.formulas.imp))
    names(pList) <- paste(ph.formulas[i],p.formulas.trans[j],l.formulas.imp)
    for (k in 1:length(l.formulas.imp)) {
      pList[[k]] <- fit.all.trans.models(
        l.formulas.imp[k],ph.formulas[i], p.formulas.trans[j]
        )
    }
    phiList[[j]] <- pList
  }
  outList.imp.trans[[i]] <- phiList
}

modelList.imp.trans <- unlist(
  unlist(outList.imp.trans, recursive = FALSE), recursive = FALSE
  )





# Fit the models for transects/demographics

outList.dem.trans <- vector('list', length = length(ph.formulas))

for (i in 1:length(ph.formulas)) {
  phiList <- vector('list', length = length(p.formulas.trans))
  for (j in 1:length(p.formulas.trans)) {
    pList <- vector('list', length = length(l.formulas.dem))
    names(pList) <- paste(ph.formulas[i],p.formulas.trans[j],l.formulas.dem)
    for (k in 1:length(l.formulas.dem)) {
      pList[[k]] <- fit.all.trans.models(
        l.formulas.dem[k],ph.formulas[i], p.formulas.trans[j]
        )
    }
    phiList[[j]] <- pList
  }
  outList.dem.trans[[i]] <- phiList
}

modelList.dem.trans <- unlist(
  unlist(outList.dem.trans, recursive = FALSE), recursive = FALSE
  )





# Fit the models for cameras/impervious

outList.imp.cam <- vector('list', length = length(p.formulas.cam))

for (i in 1:length(p.formulas.cam)) {
  p.list <- vector('list', length = length(l.formulas.imp))
  names(p.list) <- paste(p.formulas.cam[i], l.formulas.imp)
  for (j in 1:length(l.formulas.imp)) {
    p.list[[j]] <- fit.all.cam.models(p.formulas.cam[i], l.formulas.imp[j])
  }
  outList.imp.cam[[i]] <- p.list
}

modelList.imp.cam <- unlist(outList.imp.cam, recursive = FALSE)





# Fit the models for cameras/demographics

outList.dem.cam <- vector('list', length = length(p.formulas.cam))

for (i in 1:length(p.formulas.cam)) {
  p.list <- vector('list', length = length(l.formulas.dem))
  names(p.list) <- paste(p.formulas.cam[i], l.formulas.dem)
  for (j in 1:length(l.formulas.dem)) {
    p.list[[j]] <- fit.all.cam.models(p.formulas.cam[i], l.formulas.dem[j])
  }
  outList.dem.cam[[i]] <- p.list
}

modelList.dem.cam <- unlist(outList.dem.cam, recursive = FALSE)





# Make AICc tables

AICtab.imp.trans <- aictab(
  cand.set = modelList.imp.trans, modnames = names(modelList.imp.trans)
)

AICtab.dem.trans <- aictab(
  cand.set = modelList.dem.trans, modnames = names(modelList.dem.trans)
)

AICtab.imp.cam <- aictab(
  cand.set = modelList.imp.cam, modnames = names(modelList.imp.cam)
)

AICtab.dem.cam <- aictab(
  cand.set = modelList.dem.cam, modnames = names(modelList.dem.cam)
)



# ---------------------------------------------------------------------------------*
# ---- Format model results ----
# ---------------------------------------------------------------------------------*

# Function to weight predicted abundance at each site by model AICc weights

weightPredictions <- function(modResults, AICTab, method) {
  
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

unweightAbund.imp.trans <- lapply(
  modelList.imp.trans, predict, type = 'lambda', appenddata = TRUE
)

unweightAbund.dem.trans <- lapply(
  modelList.dem.trans, predict, type = 'lambda', appenddata = TRUE
)

unweightAbund.imp.cam <- lapply(
  modelList.imp.cam, predict, type = 'state', appenddata = TRUE
)

unweightAbund.dem.cam <- lapply(
  modelList.dem.cam, predict, type = 'state', appenddata = TRUE
)



# Weight predicted abundance by AIcC weights

weightAbund.imp.trans <- weightPredictions(
  unweightAbund.imp.trans, AICtab.imp.trans, 'transect'
)

weightAbund.dem.trans <- weightPredictions(
  unweightAbund.dem.trans, AICtab.dem.trans, 'transect'
)

weightAbund.imp.cam <- weightPredictions(
  unweightAbund.imp.cam, AICtab.imp.cam, 'camera'
)

weightAbund.dem.cam <- weightPredictions(
  unweightAbund.dem.cam, AICtab.dem.cam, 'camera'
)


# Write dataframes to CSV

write.csv(weightAbund.imp.trans, 'data/abundImpTrans.csv', row.names = FALSE)
write.csv(weightAbund.dem.trans, 'data/abundDemTrans.csv', row.names = FALSE)
write.csv(weightAbund.imp.cam, 'data/abundImpCam.csv', row.names = FALSE)
write.csv(weightAbund.dem.cam, 'data/abundDemCam.csv', row.names = FALSE)



