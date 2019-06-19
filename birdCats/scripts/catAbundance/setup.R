# functions ---------------------------------------------------------------

# Function searches packages in installed package list, installs them if they are not present, and loads the library:

smartLibrary <-
  function(packageVector) {
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


# Function that scales variables:

scaleVar <-
  function(var) {
    (var - mean(var, na.rm = TRUE)) / sd(var, na.rm = TRUE)
  }

# Function that scales and transposes detection covariates to wide format:

transposeCovariate <-
  function(data, covariate) {
    transMat <- data %>%
      select(site, visit, cov = covariate) %>%
      mutate(cov = scaleVar(cov)) %>%
      spread(visit, cov) %>%
      select(-site) %>%
      as.matrix()
    colnames(transMat) <- NULL
    return(transMat)
  }

# load libraries ----------------------------------------------------------

smartLibrary(
  c('unmarked', 'AICcmodavg','MuMIn', 'tidyverse')
)


# load data ---------------------------------------------------------------

options(stringsAsFactors = F)


# Visit data (time and weather):

catSiteActivity <- 
    read_csv('data/catDataActivity.csv')


# Camera sampling data:

catCam <- 
  read_csv('data/catDataCamera.csv') %>%
  filter(!is.na(species), !is.na(cameraID))


# Total Neighborhood Nestwatch site list:

catSites <- 
  read_csv('data/catSiteData.csv') %>%
  select(site)


# Transect sampling data:

catTransect <- 
  read_csv('data/catDataTransect.csv') %>%
  filter(
    !is.na(count),
    species == 'cat'
    ) %>%
  left_join(
    catSiteActivity %>%
      filter(activity == 'transect'),
    by = c('site', 'visit')
  ) %>%
  mutate(
    visit = as.factor(visit),
    doy = as.numeric(strftime(date, format = '%j'))) %>%
  arrange(site)

# Education data by site

eduData <-
  read_csv('data/eduHS.csv')


# Covariate data by site

covs <- 
  read_csv('data/covariateData.csv') %>%
  select(-c(eduHS, eduC)) %>%
  left_join(eduData) %>%
  tbl_df %>%
  arrange(site)


# Camera detection history of cats for each site:

umfCam <- 
  read_csv('data/catCamDetection.csv')


# Scaled camera detection covariates

camDetCovs <- 
  read_csv('data/camDetCovs.csv') %>%
  mutate(doy = date %>%
           strftime(format = '%j') %>%
           as.numeric) %>%
  mutate_at(
    c('tempHigh', 'tempLow', 'dewLow', 'dewHigh', 'doy'),
    function(x){scaleVar(x)}
  )
