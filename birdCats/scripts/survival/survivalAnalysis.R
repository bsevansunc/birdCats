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
  filter(
    age %in% c('AHY', 'SY', 'ASY'),
    sex %in% c('M','F'),
    !is.na(bci)) %>%
  as.data.frame() %>%
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
# Issue with this that I'll fix later--if you scale inside the full
# frame, I think the means and stdevs will be off.

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
                          'can', 'cat','imp2')),
    p = list(static = c('sex')))
  # Data input:
  proc <- getSppFrame(data, spp, catVar, scaled) %>%
    data.frame %>%
    process.data(groups="sex")
  ddl <- make.design.data(proc, design.parameters)
  
  # Models for Phi parameter:
  Phi.dot <- list(formula = ~1)
  # Additive, no cats:
  Phi.1 <- list(formula = ~sex + bci)
  Phi.2 <- list(formula = ~sex + bci + imp)
  Phi.3 <- list(formula = ~sex + bci + imp + imp2)
  Phi.4 <- list(formula = ~sex + bci + can)
  Phi.5 <- list(formula = ~sex + bci + imp + can)
  Phi.6 <- list(formula = ~sex + bci + imp + can + imp2)
  # Cat additive
  Phi.7 <- list(formula = ~sex + bci + cat)
  Phi.8 <- list(formula = ~sex + bci + cat + imp)
  Phi.9 <- list(formula = ~sex + bci + cat + imp + imp2)
  Phi.10 <- list(formula = ~sex + bci + cat + can)
  Phi.11 <- list(formula = ~sex + bci + cat + imp + can)
  Phi.12 <- list(formula = ~sex + bci + cat + imp + can + imp2)
  # sex cat interaction
  Phi.13 <- list(formula = ~sex*cat + bci)
  Phi.14 <- list(formula = ~sex*cat + bci + imp)
  Phi.15 <- list(formula = ~sex*cat + bci + imp + imp2)
  Phi.16 <- list(formula = ~sex*cat + bci + can)
  Phi.17 <- list(formula = ~sex*cat + bci + imp + can)
  Phi.18 <- list(formula = ~sex*cat + bci + imp + can + imp2)
  # bci cat interaction
  Phi.19 <- list(formula = ~bci*cat + sex)
  Phi.20 <- list(formula = ~bci*cat + sex + imp)
  Phi.21 <- list(formula = ~bci*cat + sex + imp + imp2)
  Phi.22 <- list(formula = ~bci*cat + sex + can)
  Phi.23 <- list(formula = ~bci*cat + sex + imp + can)
  Phi.24 <- list(formula = ~bci*cat + sex + imp + can + imp2)
  # imp cat interaction
  Phi.25 <- list(formula = ~imp*cat + sex + bci)
  Phi.26 <- list(formula = ~imp*cat + sex + bci + imp2)
  Phi.27 <- list(formula = ~imp*cat + sex + bci + can)
  Phi.28 <- list(formula = ~imp*cat + sex + bci + can + imp2)
  # can cat interaction
  Phi.29 <- list(formula = ~can*cat + sex + bci)
  Phi.30 <- list(formula = ~can*cat + sex + bci + imp)
  Phi.31 <- list(formula = ~can*cat + sex + bci + imp + imp2)
  # imp*can interactions
  Phi.32 <- list(formula = ~imp*can + sex + bci)
  Phi.33 <- list(formula = ~imp*can + sex + bci + imp2)
  # imp*can interactions plus cats
  Phi.34 <- list(formula = ~imp*can + sex + bci + cat)
  Phi.35 <- list(formula = ~imp*can + sex + bci + cat + imp2)
  Phi.36 <- list(formula = ~imp*can*cat + sex + bci)
  Phi.37 <- list(formula = ~imp*can*cat + sex + bci + imp2)
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
    filter(data, !site %in% removeSites),
    spp = sppVector[i],
    catVar = 'cCats',
    scaled = FALSE)
  modelListSppTransect[[i]] <- fitModels(
    data,
    spp = sppVector[i],
    catVar = 'tCats',
    scaled = FALSE)
}


# Function that calculates the weighted betas for a given sp and cov

weightedBeta <- function(modellist, spp, cov){
  betas <- vector('numeric',length=76)
  for(i in 1:76){
    betas[i] <- modellist[[spp]][[i]]$results$beta$Phi[cov]
  }
  
  table <- model.table(modellist[[spp]])
  table$index <- as.numeric(row.names(table))
  table <- table[order(table$index),]
  
  return(sum(na.omit(table$weight*betas)[1:length(na.omit(betas))]))
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


weightedBeta(modelListSppCamera, 'SOSP', 'cat')






