source('scripts/catAbundance/setup.R')


# prepare data for camera analysis --------------------------------------

# Single out transect-only sites:

removeSites <-
  c('OLONMARDC1',
    'WOLFKARDC1',
    'WOLFAMYDC1',
    'GERYERIMD1',
    'MISSEDDC1')

# Create a dataframe of potential cat detection covariates by site:

camCovs <- 
  covs %>%
  select(site) %>%
  bind_cols(
    covs %>%
      select(can:eduHS, -marred) %>%
      mutate_all(scaleVar)) %>%
  filter(!site %in% removeSites)

# Create an unmarkedFramePCount object for pcount

camUmfWithCovs <- 
  unmarkedFramePCount(
    umfCam[, -1],
    siteCovs = camCovs,
    obsCovs = camDetCovs
  )


# camera null model fitting ---------------------------------------------

# Detection formulas:

p_formulas <-
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
    '~tempHigh + tempLow',
    '~1'
  )

# Fit the models:

p_mods <-
  map(
    p_formulas,
    function(x){
      pcount(
        formula = as.formula(paste(x, '~1')),
        data = camUmfWithCovs,
        mixture = 'NB',
        K = 50)
    }) %>%
  set_names(p_formulas)

aictab(p_mods)

# impervious surface -------------------------------------------------

imp_formulas <-
  c('~imp', '~imp + I(imp^2)', '~1')

imp_mods <-
  map(
    imp_formulas,
    function(x){
      pcount(
        formula = as.formula(paste('~dewHigh', x)),
        data = camUmfWithCovs,
        mixture = 'NB',
        K = 50)
    }) %>%
  set_names(imp_formulas)

aictab(imp_mods)




# human demography ------------------------------------------------

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
      pcount(
        formula = as.formula(paste('~dewHigh', x)),
        data = camUmfWithCovs,
        mixture = 'NB',
        K = 50)
    }) %>%
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
      pcount(
        formula = as.formula(paste('~dewHigh', x)),
        data = camUmfWithCovs,
        mixture = 'NB',
        K = 50)
    }) %>%
  set_names(hDem_formulas_reduced)

aictab(hDem_mods_reduced)
