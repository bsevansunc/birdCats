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
  select(site) %>%
  bind_cols(
    covs %>%
      select(can:eduHS) %>%
      mutate_all(scaleVar)
  )

# Create a dataframe of potential cat detection covariates by site:

transDetCovs <- 
  catTransect %>%
  select(-c(species:date)) %>%
  distinct %>%
  mutate_at(
    c('time', 'temp', 'dew', 'doy'),
    scaleVar)

# Create a list of wide-form detection covariates:

yCov_names <- 
  c('time', 'temp','dew', 'doy')

ySiteCovs <-
  map(
    yCov_names,
    function(x){
      transposeCovariate(transDetCovs, x)
    }
  ) %>%
  set_names(yCov_names)

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

# distance function -------------------------------------------------------

# Following vignette: https://rstudio-pubs-static.s3.amazonaws.com/221408_23c61679859e48e6bae0b9c5c2e48a92.html

# Determine which distance function best describes how counts vary by distance:

detection_functions <-
  c('halfnorm', 'hazard', 'exp', 'uniform')

distance_mods <-
  map(
    detection_functions,
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = '~1',
        pformula = '~1',
        data = gUmfWithCovs,
        keyfun = x,
        mixture = 'NB',
        K = 50,
        output = 'abund'
      )
    }
  ) %>% 
  set_names(detection_functions) 

aictab(distance_mods)

# Hazard distance function was best supported.

# availability and detection ----------------------------------------------

# formulas for phi and p parameters:

phiP_formulas <-
  crossing(
    phi = c('~1', '~time', '~I(time^2)'),
    p =  c(
      '~temp',
      '~temp + time',
      '~temp + dew',
      '~temp + time + dew',
      '~temp + time + I(time^2) + dew',
      '~temp*dew + time',
      '~temp*dew + time + I(time^2)',
      '~time',
      '~time + I(time^2)',
      '~time + dew',
      '~time + I(time^2) + dew',
      '~dew',
      '~1'
    )
  )

distance_mods <-
  map(
    1:nrow(phiP_formulas),
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = phiP_formulas$phi[x],
        pformula = phiP_formulas$p[x],
        data = gUmfWithCovs,
        keyfun = 'hazard',
        mixture = 'NB',
        K = 50,
        output = 'abund'
      )
    }
  ) %>% 
  set_names(
    paste(
      phiP_formulas$phi, 
      phiP_formulas$p))

aictab(distance_mods)

# impervious surface ------------------------------------------------------

imp_formulas <-
  c('~imp', '~imp + I(imp^2)', '~1')

imp_mods <-
  map(
  imp_formulas,
  function(x){
    gdistsamp(
      lambdaformula = x,
      phiformula = '~time',
      pformula = '~1',
      data = gUmfWithCovs,
      keyfun = 'hazard',
      mixture = 'NB',
      K = 50,
      output = 'abund'
    )
  }
) %>% 
  set_names(imp_formulas)


aictab(imp_mods)

# Conduct a goodness-of-fit test

Nmix.gof.test(imp_mods[[1]])
