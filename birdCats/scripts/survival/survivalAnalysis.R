library('tidyverse')
library('marked')
library('stringr')

removeSites <- c('GERYERIMD1', 'OLONMARDC1', 'MISSEDDC1', 'WOLFKARDC1', 'WOLFAMYDC1')

# Read in data

data <- read_csv('data/ch.csv')


# Filter out individuals of unknown sex and body condition

bird.c <- data %>%
  filter(
    !is.na(bci),
    age != "U" & age != 'HY',
    sex != 'U',
    !site %in% removeSites) %>%
  as.data.frame


# Scales covs; for testing scaled versus unscaled
# (messy, I know, but they do have to be scaled outside the full dataframe)

cCatsScale <- bird.c %>% select(site,cCats) %>% unique
cCatsScale$cCats <- scale(cCatsScale$cCats)[,1]
bird.c <- bird.c %>% select(-cCats)
bird.c <- left_join(bird.c,cCatsScale,by='site')

tCatsScale <- bird.c %>% select(site,tCats) %>% unique
tCatsScale$tCats <- scale(tCatsScale$tCats)[,1]
bird.c <- bird.c %>% select(-tCats)
bird.c <- left_join(bird.c,tCatsScale,by='site')

impScale <- bird.c %>% select(site,imp) %>% unique
impScale$imp <- scale(impScale$imp)[,1]
bird.c <- bird.c %>% select(-imp)
bird.c <- left_join(bird.c,impScale,by='site')

imp2Scale <- bird.c %>% select(site,imp2) %>% unique
imp2Scale$imp2 <- scale(imp2Scale$imp2)[,1]
bird.c <- bird.c %>% select(-imp2)
bird.c <- left_join(bird.c,imp2Scale,by='site')



# Set sex as a factor with levels 'M' and 'F'

bird.c$sex <- as.factor(bird.c$sex)


# By species

grca.c <- bird.c %>% filter(species == 'GRCA')


# Process the data

grca.c.proc <- process.data(grca.c, groups = "sex")


# Design parameters
design.Phi <- list(static=c('imp','sex','bci','imp2', 'can', 'can2', 'tCats','cCats'))
design.p <- list(static = c('sex'))
design.parameters <- list(Phi = design.Phi, p = design.p)

grca.c.ddl <- make.design.data(grca.c.proc, parameters = design.parameters)




fit.models.c <- function(proc,ddl){
  Phi.dot <- list(formula = ~1)
  # Additive
  Phi.1 <- list(formula = ~sex)
  Phi.2 <- list(formula = ~sex + bci)
  Phi.3 <- list(formula = ~sex + bci + imp)
  Phi.4 <- list(formula = ~sex + bci + imp + imp2)
  # Cat additive
  Phi.5 <- list(formula = ~sex + cCats)
  Phi.6 <- list(formula = ~sex + bci + cCats)
  Phi.7 <- list(formula = ~sex + bci + cCats + imp)
  Phi.8 <- list(formula = ~sex + bci + cCats + imp +imp2)
  
  # Interactions
  Phi.9 <- list(formula = ~sex*bci)
  Phi.10 <- list(formula = ~sex*bci + imp)
  Phi.11 <- list(formula = ~sex*bci + imp + imp2)
  Phi.12 <- list(formula = ~bci + imp*sex)
  Phi.13 <- list(formula = ~bci + imp*sex + imp2)
  Phi.14 <- list(formula = ~sex*bci + imp*sex)
  Phi.15 <- list(formula = ~sex*bci + imp*sex + imp2)
  # Cat interactions
  Phi.16 <- list(formula = ~sex*bci + cCats)
  Phi.17 <- list(formula = ~sex*bci + imp*cCats)
  Phi.18 <- list(formula = ~sex*bci + imp*cCats + imp2)
  Phi.19 <- list(formula = ~bci + sex*cCats)
  Phi.20 <- list(formula = ~sex + bci*cCats)
  Phi.21 <- list(formula = ~sex*cCats + bci*cCats)
  Phi.22 <- list(formula = ~bci + sex*cCats + imp*cCats)
  Phi.23 <- list(formula = ~sex + bci*cCats + imp*cCats)
  Phi.24 <- list(formula = ~sex*cCats + bci*cCats + imp*cCats)
  Phi.25 <- list(formula = ~bci + sex*cCats + imp*cCats + imp2)
  Phi.26 <- list(formula = ~sex + bci*cCats + imp*cCats + imp2)
  Phi.27 <- list(formula = ~sex*cCats + bci*cCats + imp*cCats + imp2)
  Phi.28 <- list(formula = ~bci + sex*cCats + imp*sex + cCats*imp)
  Phi.29 <- list(formula = ~sex + bci*cCats + imp*sex + cCats*imp)
  Phi.30 <- list(formula = ~sex*cCats + bci*cCats + imp*sex + cCats*imp)
  Phi.31 <- list(formula = ~bci + sex*cCats + imp*sex + cCats*imp + imp2)
  Phi.32 <- list(formula = ~sex + bci*cCats + imp*sex + cCats*imp + imp2)
  Phi.33 <- list(formula = ~sex*cCats + bci*cCats + imp*sex + cCats*imp + imp2)
  Phi.34 <- list(formula = ~sex*bci + sex*cCats + bci*cCats + imp*cCats)
  Phi.35 <- list(formula = ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp2)
  Phi.36 <- list(formula = ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp*sex)
  Phi.37 <- list(formula = ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp*sex*imp2)
  
  p.sex <- list(formula = ~sex)
  
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml,data=proc,ddl=ddl,external=FALSE,accumulate=FALSE)
  return(results)
}


# Fit the models

grca.models.c <- fit.models.c(grca.c.proc,grca.c.ddl)




