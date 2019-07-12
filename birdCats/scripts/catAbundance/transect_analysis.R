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

yCov_names <- 
  c('time', 'temp','dew', 'doy')

transDetCovs <- 
  catTransect %>%
  select(-c(species:date)) %>%
  distinct %>%
  mutate_at(
    yCov_names,
    scaleVar)

# Create a list of wide-form detection covariates:

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

# probability density function -------------------------------------------------------

# Following vignette: https://rstudio-pubs-static.s3.amazonaws.com/221408_23c61679859e48e6bae0b9c5c2e48a92.html

# Determine which distance function best describes how counts vary by distance:

detection_functions <-
  c('halfnorm', 'hazard')

phi_formulas <- '~1'
p_formulas <- '~1'

formula_frame <-
  crossing(detection_functions, phi_formulas, p_formulas)

detection_mods <-
  map(
    1:nrow(formula_frame),
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = formula_frame[x,]$phi_formulas,
        pformula = formula_frame[x,]$p_formulas,
        data = gUmfWithCovs,
        keyfun = formula_frame[x,]$detection_functions,
        mixture = 'NB',
        K = 50,
        output = 'abund')
    }) %>%
  set_names(detection_functions)

aictab(detection_mods)

# Hazard distance function was best supported.

# availability and detection functions ------------------------------------

detection_functions <- 'hazard'

phi_formulas <-
  c('~1',
    '~time',
    '~doy',
    '~time + doy',
    '~time + I(time^2)',
    '~doy + I(doy^2)'
  )

p_formulas <-
  c(
    '~1',
    '~temp',
    '~time',
    '~dew',
    '~temp + time',
    '~temp + dew',
    '~time + dew',
    '~temp + time + dew',
    '~time + I(time^2)',
    '~temp + time + I(time^2)',
    '~time + I(time^2) + dew',
    '~temp + I(time^2) + time + dew'
  )

formula_frame <-
  crossing(detection_functions, phi_formulas, p_formulas)

phi_p_mods <-
  map(
    1:nrow(formula_frame),
    function(x){
      gdistsamp(
        lambdaformula = '~ 1',
        phiformula = formula_frame[x,]$phi_formulas,
        pformula = formula_frame[x,]$p_formulas,
        data = gUmfWithCovs,
        keyfun = formula_frame[x,]$detection_functions,
        mixture = 'NB',
        K = 50,
        output = 'abund')
      }) %>%
  set_names(
    str_c(
      formula_frame$detection_functions,
      formula_frame$phi_formulas,
      formula_frame$p_formulas,
      sep = ' '
    ))


aictab(phi_p_mods)

# evaluate overdispersion -------------------------------------------------

# Create a global model (most complex ... for abundance parameter): 

global_mod <-
  gdistsamp(
    lambdaformula = '~ imp + I(imp^2)',
    phiformula = '~doy + I(doy^2)',
    pformula = '~1',
    data = gUmfWithCovs,
    keyfun = 'halfnorm',
    mixture = 'NB',
    K = 50,
    output = 'abund')

# Evaluate overdispersion (c-hat):

Nmix.gof.test(global_mod)

# availability ------------------------------------------------------------

phi_formulas <-
  c('~1',
    # 1 variable:
    '~temp',
    '~time',
    '~dew',
    '~doy',
    # 2 variables:
    '~temp + time',
    '~temp + dew',
    '~temp + doy',
    '~time + dew',
    '~time + doy',
    '~dew + doy',
    # 3 variables:
    '~temp + time + dew',
    '~temp + time + doy',
    '~temp + dew + doy',
    '~time + dew + doy',
    # 4 variables:
    '~temp + time + dew + doy'
  )

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

# detectability -----------------------------------------------------------

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

# availability and detectability ------------------------------------------

phi_formulas <-
  c('~1', '~time')

p_formulas <-
  c(
    '~temp',
    '~temp + time',
    '~temp + dew',
    '~temp + time + dew',
    # '~temp + doy',
    # '~temp + time + dew + doy',
    '~time',
    '~time + dew',
    # '~time + dew + doy',
    '~dew',
    # '~dew + doy',
    # '~doy',
    '~1'
  )

formula_frame <-
  crossing(phi_formulas, p_formulas)

phiP_mods <-
  map(
    1:nrow(formula_frame),
    function(x){
        gdistsamp(
          lambdaformula = '~ 1',
          phiformula = formula_frame[x, 'phi_formulas'],
          pformula = formula_frame[x, 'p_formulas'],
          data = gUmfWithCovs,
          keyfun = 'hazard',
          mixture = 'NB',
          K = 50,
          output = 'abund')
        }) %>% 
  set_names(
    str_c(
      formula_frame$phi_formulas,
      formula_frame$p_formulas,
      sep = ' '
  ))

aictab(phiP_mods)

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

# human demography --------------------------------------------------------

hDem_formulas <-
  c(
    '~1',
    # One variable:
    '~medianIncome',
    '~hDensity',
    '~age',
    '~marred',
    '~eduHS',
    # Two variables:
    '~medianIncome + hDensity',
    '~medianIncome + age',
    '~medianIncome + marred',
    '~medianIncome + eduHS',
    '~hDensity + age',
    '~hDensity + marred',
    '~hDensity + eduHS',
    '~age + marred',
    '~age + eduHS',
    '~marred + eduHS',
    # Three variables:
    '~medianIncome + hDensity + age',
    '~medianIncome + hDensity + marred',
    '~medianIncome + hDensity + eduHS',
    '~medianIncome + age + marred',
    '~medianIncome + age + eduHS',
    '~medianIncome + marred + eduHS',
    '~hDensity + age + marred',
    '~hDensity + age + eduHS',
    '~hDensity + marred + eduHS',
    '~age + marred + eduHS',
    # Four variables:
    '~medianIncome + hDensity + age + marred',
    '~medianIncome + hDensity + age + eduHS',
    '~hDensity + age + marred + eduHS',
    # Five variables:
    '~medianIncome + hDensity + age + marred + eduHS',
    # With quadratic hDensity, 1 variable:
    '~hDensity + I(hDensity^2)',
    # With quadratic hDensity, 2 variables:
    '~medianIncome + hDensity + I(hDensity^2)',
    '~hDensity + I(hDensity^2) + age',
    '~hDensity + I(hDensity^2) + marred',
    '~hDensity + I(hDensity^2) + eduHS',
    # With quadratic hDensity, 3 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age',
    '~medianIncome + hDensity + I(hDensity^2) + marred',
    '~medianIncome + hDensity + I(hDensity^2) + eduHS',
    '~hDensity + I(hDensity^2) + age + marred',
    '~hDensity + I(hDensity^2) + age + eduHS',
    '~hDensity + I(hDensity^2) + marred + eduHS',
    # With quadratic hDensity, 4 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age + marred',
    '~medianIncome + hDensity + I(hDensity^2) + age + eduHS',
    '~hDensity + I(hDensity^2) + age + marred + eduHS',
    # With quadratic hDensity, 5 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age + marred + eduHS',
    # With quadratic income, 1 variable:
    '~medianIncome + I(medianIncome^2)',
    # With quadratic income, 2 variables:
    '~medianIncome + I(medianIncome^2) + hDensity',
    '~medianIncome + I(medianIncome^2) + age',
    '~medianIncome + I(medianIncome^2) + marred',
    '~medianIncome + I(medianIncome^2) + eduHS',
    # With quadratic income, 3 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + age',
    '~medianIncome + I(medianIncome^2) + hDensity + marred',
    '~medianIncome + I(medianIncome^2) + hDensity + eduHS',
    '~medianIncome + I(medianIncome^2) + age + marred',
    '~medianIncome + I(medianIncome^2) + age + eduHS',
    '~medianIncome + I(medianIncome^2) + marred + eduHS',
    # With quadratic income, 4 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + age + marred',
    '~medianIncome + I(medianIncome^2) + hDensity + age + eduHS',
    # With quadratic income, 5 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + age + marred + eduHS',
    # Quadratic income and hDensity, 2 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2)',
    # Quadratic income and hDensity, 3 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + age',
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + marred',
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + eduHS',
    # Quadratic income and hDensity, 4 variables:
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + age + marred',
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + age + eduHS',
    # Quadratic income and hDensity, 5 variables (global mod):
    '~medianIncome + I(medianIncome^2) + hDensity + I(hDensity^2) + age + marred + eduHS'
    )

hDem_mods <-
  map(
    hDem_formulas,
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
  set_names(hDem_formulas)

aictab(hDem_mods)

# Compare models with quadratic housing density term vs. linear only:

data.frame(aictab(hDem_mods)) %>% 
  as_tibble() %>%
  mutate(
    mods = case_when(
      str_detect(Modnames, 'hDensity') &
        !str_detect(Modnames, 'hDensity\\^2') ~ 'hDensity',
      str_detect(Modnames, 'hDensity\\^2') ~ 'hDensity2',
    TRUE ~ 'other'
    )) %>%
  group_by(mods) %>%
  summarize(cumWt = sum(AICcWt))

# Compare models with quadratic income term vs. linear only:

data.frame(aictab(hDem_mods)) %>% 
  as_tibble() %>%
  mutate(
    mods = case_when(
      str_detect(Modnames, 'medianIncome') &
        !str_detect(Modnames, 'medianIncome\\^2') ~ 'medianIncome',
      str_detect(Modnames, 'medianIncome\\^2') ~ 'medianIncome2',
      TRUE ~ 'other'
    )) %>%
  group_by(mods) %>%
  summarize(cumWt = sum(AICcWt))


# human demographics, reduced ---------------------------------------------

hDem_formulas_reduced <-
  c(
    '~1',
    # One variable:
    '~medianIncome',
    '~age',
    '~marred',
    '~eduHS',
    # Two variables:
    '~medianIncome + age',
    '~medianIncome + marred',
    '~medianIncome + eduHS',
    '~age + marred',
    '~age + eduHS',
    '~marred + eduHS',
    # Three variables:
    '~medianIncome + age + marred',
    '~medianIncome + age + eduHS',
    '~medianIncome + marred + eduHS',
    '~age + marred + eduHS',
    # With quadratic hDensity, 1 variable:
    '~hDensity + I(hDensity^2)',
    # With quadratic hDensity, 2 variables:
    '~medianIncome + hDensity + I(hDensity^2)',
    '~hDensity + I(hDensity^2) + age',
    '~hDensity + I(hDensity^2) + marred',
    '~hDensity + I(hDensity^2) + eduHS',
    # With quadratic hDensity, 3 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age',
    '~medianIncome + hDensity + I(hDensity^2) + marred',
    '~medianIncome + hDensity + I(hDensity^2) + eduHS',
    '~hDensity + I(hDensity^2) + age + marred',
    '~hDensity + I(hDensity^2) + age + eduHS',
    '~hDensity + I(hDensity^2) + marred + eduHS',
    # With quadratic hDensity, 4 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age + marred',
    '~medianIncome + hDensity + I(hDensity^2) + age + eduHS',
    '~hDensity + I(hDensity^2) + age + marred + eduHS',
    # With quadratic hDensity, 5 variables:
    '~medianIncome + hDensity + I(hDensity^2) + age + marred + eduHS'
  )

hDem_mods_reduced <-
  map(
    hDem_formulas_reduced,
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
  set_names(hDem_formulas_reduced)

    
aictab(hDem_mods_reduced)
  
  




