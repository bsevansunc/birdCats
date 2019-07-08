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
        output = 'abund')}) %>% 
  set_names(detection_functions) 

aictab(distance_mods)

# Hazard distance function was best supported.

# Starting values for further modeling:

start_values <- 
  c(1.75, -2.62, 2.84,-.446)

# evaluate overdispersion -------------------------------------------------

# Create a global model (most complex): 

global_mod <-
  gdistsamp(
    lambdaformula = '~ imp + I(imp^2)',
    phiformula = '~time + doy',
    pformula = '~temp + dew + time + doy',
    data = gUmfWithCovs,
    keyfun = 'hazard',
    mixture = 'NB',
    K = 50,
    output = 'abund',
    se = FALSE)

# Evaluate overdispersion (c-hat):

Nmix.gof.test(global_mod)


# availability ------------------------------------------------------------

phi_formulas <-
  c('~1', '~time')

phi_mods <-
  map(
    phi_formulas,
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = x,
        pformula = '~1',
        data = gUmfWithCovs,
        keyfun = 'hazard',
        mixture = 'NB',
        K = 50,
        output = 'abund')}) %>% 
  set_names(phi_formulas)

aictab(phi_mods)

# detection ---------------------------------------------------------------

p_formulas <-
  c(
    '~temp',
    '~temp + time',
    '~temp + dew',
    '~temp + doy',
    '~temp + time + dew + doy',
    '~time',
    '~time + dew',
    '~time + dew + doy',
    '~dew',
    '~dew + doy',
    '~doy',
    '~1'
  )


p_mods <-
  map(
    p_formulas,
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = '~time',
        pformula = x,
        data = gUmfWithCovs,
        keyfun = 'hazard',
        mixture = 'NB',
        K = 50,
        output = 'abund')}) %>% 
  set_names(p_formulas)

aictab(p_mods)

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
      output = 'abund')}) %>% 
  set_names(imp_formulas)

aictab(imp_mods)


