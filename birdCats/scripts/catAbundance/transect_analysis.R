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

# fit availability and detection models -----------------------------------

# Following vignette: https://rstudio-pubs-static.s3.amazonaws.com/221408_23c61679859e48e6bae0b9c5c2e48a92.html

# Determine which distance function best describes how counts vary by distance:

detection_functions <-
  c('halfnorm', 'hazard', 'exp', 'uniform')

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
  set_names(detection_functions) %>%
  aictab()

# Hazard distance function was best supported.

# Determine whether availability (phi) varies by time of day and impervious:

availability_functions <-
  c('~1', '~time', '~I(time^2)')

map(
  availability_functions,
  function(x){
    gdistsamp(
      lambdaformula = '~ 1',
      phiformula = x,
      pformula = '~1',
      data = gUmfWithCovs,
      keyfun = 'hazard',
      mixture = 'NB',
      K = 50,
      output = 'abund'
    )
  }
) %>% 
  set_names(availability_functions) %>%
  aictab()

# Determine whether detection varies by date, time, or weather:

detection_functions <-
  c(
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

map(
  detection_functions,
  function(x){
    gdistsamp(
      lambdaformula = '~ 1',
      phiformula = '~1',
      pformula = x,
      data = gUmfWithCovs,
      keyfun = 'hazard',
      mixture = 'NB',
      K = 50,
      output = 'abund'
    )
  }
) %>% 
  set_names(detection_functions) %>%
  aictab()





# Formulas for availability (phi) and detection(p):

formulas <-
  crossing(
    phi = c('~time', '~time + I(time^2)', '~1'),
    p = c(
      '~doy + temp',
      '~doy + time',
      '~doy + dew',
      '~doy + time + I(time^2)',
      '~dew + time',
      '~dew + temp',
      '~dew + time + I(time^2)',
      '~time + temp',
      '~time + I(time^2) + temp',
      '~doy',
      '~time',
      '~time + I(time^2)',
      '~dew',
      '~temp',
      '~1'
    )
  ) %>%
  arrange(phi, p)

# Run models:

tModels <-
  purrr::map(
    1:nrow(formulas),
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = formulas$phi[x],
        pformula = formulas$p[x],
        data = gUmfWithCovs,
        keyfun = 'halfnorm',
        mixture = 'NB',
        K = 50,
        output = 'abund'
      )
    }
  ) %>%
  set_names(paste(formulas$phi, formulas$p))

# AIC table to compare availability/detectability models:

aictab(
  cand.set = tModels,
  modnames = names(tModels))
