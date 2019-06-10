source('scripts/catAbundance/setup.R')

# prepare data ------------------------------------------------------------

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
  select(-c(lon:fips)) # %>%
  # bind_cols(
  #   covs %>%
  #     select(can:age, eduHS) %>%
  #     mutate_all(function(x) x^2) %>%
  #     set_names('can2', 'imp2', 'age2', 'medianIncome2', 'hDensity2', 'eduHS2')
  # )
  
  
  # select(site) %>%
  # bind_cols(
  #   covs %>%
  #     select(can:eduHS) %>%
  #     mutate(
  #       can2 = can^2,
  #       imp2 = imp^2,
  #       age2 = age^2,
  #       medianIncome2 = medianIncome^2,
  #       hDensity2 = hDensity^2,
  #       eduHS2 = eduHS^2
  #     ) %>%
  #     mutate_all(scaleVar)
  # )

# Create a dataframe of potential cat detection covariates by site:

visitCovs <-
  c('time','temp','dew', 'doy')

transDetCovs <- 
  catTransect %>%
  select(-c(species:date)) %>%
  distinct %>%
  mutate(
    time2 = time^2, 
    temp2 = temp^2) %>%
  mutate_at(
    visitCovs,
    scaleVar)

# Create a list of wide-form detection covariates:

ySiteCovs <-
  map(visitCovs, function(x){
    transposeCovariate(transDetCovs, x)
  }) %>%
  set_names(visitCovs)

# Create unmarkedFrameGDS object for gdistsamp:

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
  c(
    '~time',
    # '~time + I(time^2)',
    'doy',
    '~1'
  )

# Formulas for detection:

pFormulas <-
  c(# No quadratic terms:
    '~doy',
    '~doy + temp',
    '~doy + temp + time',
    '~doy + temp + dew',
    '~doy + temp + time + dew',
    '~temp',
    '~temp + time',
    '~temp + dew',
    '~temp + time + dew',
    '~time',
    '~time + dew',
    '~dew',
    # Interaction of temp and dew:
    '~doy + temp * dew',
    '~doy + time + temp * dew',
    '~time + temp * dew',
    'temp*dew',
    # # Quadratic temp:
    # # '~doy + temp + I(temp^2)',
    # # '~doy + temp + I(temp^2) + time',
    # # '~doy + temp + I(temp^2) + dew',
    # # '~doy + temp + I(temp^2) + time + dew',
    # # '~doy + temp + I(temp^2) + time * dew',
    # # '~temp + I(temp^2)',
    # # '~temp + I(temp^2) + time',
    # # '~temp + I(temp^2) + dew',
    # # '~temp + I(temp^2) + time + dew',
    # # '~temp + I(temp^2) + time * dew',
    # # Quadratic time:
    # '~doy + temp + time + I(time^2)',
    # '~doy + temp + time + dew + I(time^2)',
    # '~doy + temp + time * dew + I(time^2)',
    # '~temp + time + I(time^2)',
    # '~temp + time + dew + I(time^2)',
    # '~temp + time * dew + I(time^2)',
    # '~time + I(time^2)',
    # '~time + dew + I(time^2)',
    # '~time * dew + I(time^2)',
    # # Quadratic temp and time:
    # # '~doy + temp + time + I(time^2) + I(temp^2)',
    # # '~doy + temp + time + dew + I(time^2) + I(temp^2)',
    # # '~doy + temp + time * dew + I(time^2) + I(temp^2)',
    # # '~temp + time + I(time^2) + I(temp^2)',
    # # '~temp + time + dew + I(time^2) + I(temp^2)',
    # # '~temp + time * dew + I(time^2) + I(temp^2)',
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

