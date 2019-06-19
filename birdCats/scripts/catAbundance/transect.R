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

# Create a dataframe of covariates by site:

transCovs <- 
  covs %>%
  select(-c(lon:fips)) %>%
  arrange(site)

# Create a dataframe of cat availability/detection covariates by site:

visitCovs <-
  c('time','temp','dew', 'doy')

transDetCovs <- 
  catTransect %>%
  arrange(site, visit) %>%
  select(site, visit, time:doy) %>%
  distinct %>%
  mutate_at(
    visitCovs,
    scaleVar)

# Create a list of wide-form detection covariates:

ySiteCovs <-
  map(visitCovs, function(x){
    transDetCovs %>% 
      select(site, visit, x) %>% 
      spread(visit, x) %>%
      select(-site) %>%
      as.matrix
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

# fit phi and p parameters ------------------------------------------------

# Formulas for availability (phi) and detection (p):

formulas <-
  crossing(
    phi = c(
      '~time',
      '~1'
    ), 
    p = c(
      '~temp',
      '~dew',
      '~temp + dew',
      '~1'
    )) %>%
  mutate(f = paste(phi, p))

# Build null abundance models with availability/detection covariates:

modelList <-
  map(
  1:nrow(formulas),
  function(x){
    gdistsamp(
      lambdaformula = '~ 1',
      phiformula = formulas[x,'phi'],
      pformula = formulas[x,'p'],
      data = gUmfWithCovs,
      keyfun = 'halfnorm',
      mixture = 'NB',
      K = 50
    )
  }) %>%
  set_names(formulas$f)

# AIC table to compare availability/detection models:

aictab(modelList)

# abundance model ---------------------------------------------------------

formulas <-
  crossing(
    phi = c(
      '~time'
    ), 
    p = c(
      '~dew'
    ),
    lambda = c(
      '~ 1',
      '~ imp + I(imp^2)',
      '~ hDensity',
      '~ medianIncome',
      '~ eduHS',
      '~ hDensity + medianIncome',
      '~ hDensity + eduHS',
      '~ hDensity + medianIncome+eduHS',
      '~ medianIncome + eduHS',
      '~ imp + I(imp^2) + hDensity',
      '~ imp + I(imp^2) + medianIncome',
      '~ imp + I(imp^2) + eduHS',
      '~ imp + I(imp^2) + hDensity+medianIncome',
      '~ imp + I(imp^2) + hDensity+eduHS',
      '~ imp + I(imp^2) + hDensity+medianIncome+eduHS',
      '~ imp + I(imp^2) + medianIncome+eduHS'
    )) %>%
  mutate(f = paste(phi, p, lambda))



modelList_abund <-
  map(
    1:nrow(formulas),
    function(x){
      gdistsamp(
        lambdaformula = formulas[x, 'lambda'],
        phiformula = formulas[x,'phi'],
        pformula = formulas[x,'p'],
        data = gUmfWithCovs,
        keyfun = 'halfnorm',
        mixture = 'NB',
        K = 50
      )
    }) %>%
  set_names(formulas$f)


# Create an AIC table

aictab(modelList_abund)


