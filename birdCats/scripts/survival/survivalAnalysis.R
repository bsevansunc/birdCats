#--------------------------------------------------------------------*
# ---- Set up ----
#--------------------------------------------------------------------*

library('tidyverse')
library('marked')
library('stringr')

removeSites <- c('GERYERIMD1', 'OLONMARDC1', 'MISSEDDC1',
                 'WOLFKARDC1', 'WOLFAMYDC1')

# Read and format data:

data <- read_csv('data/ch.csv') %>%
  tbl_df %>%
  filter(
    age %in% c('AHY', 'SY', 'ASY'),
    sex %in% c('M', 'F'),
    !site %in% removeSites)  %>%
  na.omit %>%
  mutate(sex = as.factor(sex)) %>%
  select(ch, bandYear, age, sex, species, bci:tCats)

#--------------------------------------------------------------------*
# ---- Functions ----
#--------------------------------------------------------------------*

# Function to scale covariates:

scaleCov <- function(cov){
  scale(cov)[,1]
}

# Function to get frame for a given species:

getSppFrame <- function(data, spp, catVar, scaled = TRUE){
  sppData <- data %>%
    filter(species == spp) %>%
    mutate_(cat = catVar) %>%
    select(-c(cCats, tCats))
  if(scaled == TRUE){
    sppData <- bind_cols(
      sppData %>%
        select(ch:bci),
      sppData %>%
        select(imp:tCats) %>%
        mutate_all(scaleCov)
    )
  }
  return(sppData)
}

# Function to fit models, camera:

fitModels <- function(data, spp, catVar, scaled = FALSE){
  # Design parameters:
  design.parameters <- list(
    Phi = list(static = c('imp','sex','bci',
                          'can', 'cat')),
    p = list(static = c('sex')))
  # Data input:
  proc <- getSppFrame(data, spp, catVar, scaled) %>%
    data.frame %>%
    process.data
  ddl <- make.design.data(proc, design.parameters)
  # Models for Phi parameter:
  Phi.dot <- list(formula = ~1)
  # Additive, no cats:
  Phi.1 <- list(formula = ~sex + bci)
  Phi.2 <- list(formula = ~sex + bci + imp)
  Phi.3 <- list(formula = ~sex + bci + can)
  Phi.4 <- list(formula = ~sex + bci + imp + can)
  # Cat additive
  Phi.5 <- list(formula = ~sex + bci + cat)
  Phi.5 <- list(formula = ~sex + bci + cat + imp)
  Phi.6 <- list(formula = ~sex + bci + cat + can)
  Phi.7 <- list(formula = ~sex + bci + cat + imp + can)
  # sex cat interaction
  Phi.8 <- list(formula = ~sex*cat + bci)
  Phi.9 <- list(formula = ~sex*cat + bci + imp)
  Phi.10 <- list(formula = ~sex*cat + bci + can)
  Phi.11 <- list(formula = ~sex*cat + bci + imp + can)
  # bci cat interaction
  Phi.11 <- list(formula = ~bci*cat + sex)
  Phi.12 <- list(formula = ~bci*cat + sex + imp)
  Phi.13 <- list(formula = ~bci*cat + sex + can)
  Phi.14 <- list(formula = ~bci*cat + sex + imp + can)
  # imp cat interaction
  Phi.15 <- list(formula = ~imp*cat + sex + bci)
  Phi.16 <- list(formula = ~imp*cat + sex + bci + can)
  # can cat interaction
  Phi.17 <- list(formula = ~can*cat + sex + bci)
  Phi.18 <- list(formula = ~can*cat + sex + bci + imp)
  # imp*can interactions
  Phi.19 <- list(formula = ~imp*can + sex + bci)
  Phi.20 <- list(formula = ~imp*can + sex + bci + cat)
  Phi.21 <- list(formula = ~imp*can*cat + sex + bci)
  # Model for p parameter
  p.dot <- list(formula = ~1)
  p.sex <- list(formula = ~sex)
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml,data=proc,ddl=ddl,
                         external=FALSE,accumulate=FALSE)
  return(results)
}


#--------------------------------------------------------------------*
# ---- Model running ----
#--------------------------------------------------------------------*

# Fit the models:

sppVector <- data$species %>% unique %>% sort

modelListSppCamera <- vector('list', length = length(sppVector))
modelListSppTransect <- vector('list', length = length(sppVector))

names(modelListSppCamera) <- sppVector
names(modelListSppTransect) <- sppVector


for(i in 1:length(sppVector)){
  modelListSppCamera[[i]] <- fitModels(
    data,
    spp = sppVector[i],
    catVar = 'cCats',
    scaled = FALSE)
  modelListSppTransect[[i]] <- fitModels(
    data,
    spp = sppVector[i],
    catVar = 'tCats',
    scaled = FALSE)
}

# Example output:

modelListSppCamera$GRCA

modelListSppCamera$GRCA$Phi.21.p.dot

modelListSppCamera$SOSP$Phi.7.p.sex

modelListSppCamera$SOSP$Phi.7.p.sex$results$reals$Phi %>%
  ggplot(aes(x = cat, y = estimate)) +
  geom_point(cex = 3, alpha = .5) +
  theme_bw()

modelListSppCamera$GRCA$Phi.7.p.sex$results$reals$Phi %>%
  ggplot(aes(x = cat, y = estimate)) +
  geom_point(cex = 3, alpha = .5) +
  theme_bw()











