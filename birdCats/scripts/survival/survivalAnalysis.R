library('tidyverse')
library('marked')
library('stringr')

removeSites <- c('GERYERIMD1', 'OLONMARDC1', 'MISSEDDC1', 'WOLFKARDC1', 'WOLFAMYDC1')

# Read in data

data <- read_csv('data/ch.csv')


# Filter out individuals of unknown sex and body condition
bird.t <- data %>%
  filter(
    !is.na(bci),
    age != "U" & age != 'HY',
    sex != 'U') %>%
  as.data.frame

bird.c <- data %>%
  filter(
    !is.na(bci),
    age != "U" & age != 'HY',
    sex != 'U',
    !site %in% removeSites) %>%
  as.data.frame


# # Scales cCats; for testing scaled versus unscaled:
# test <- bird.c %>% select(site,cCats) %>% unique
# test$cCats <- scale(test$cCats)[,1]
# bird.c <- bird.c %>% select(-cCats)
# bird.c <- left_join(bird.c,test,by='site')



# Set sex as a factor with levels 'M' and 'F'
bird.t$sex <- as.factor(bird.t$sex)
bird.c$sex <- as.factor(bird.c$sex)


# By species
amro.t <- bird.t %>% filter(species == 'AMRO')
cach.t <- bird.t %>% filter(species == 'CACH')
carw.t <- bird.t %>% filter(species == 'CARW')
grca.t <- bird.t %>% filter(species == 'GRCA')
howr.t <- bird.t %>% filter(species == 'HOWR')
noca.t <- bird.t %>% filter(species == 'NOCA')
sosp.t <- bird.t %>% filter(species == 'SOSP')

amro.c <- bird.c %>% filter(species == 'AMRO')
cach.c <- bird.c %>% filter(species == 'CACH')
carw.c <- bird.c %>% filter(species == 'CARW')
grca.c <- bird.c %>% filter(species == 'GRCA')
howr.c <- bird.c %>% filter(species == 'HOWR')
noca.c <- bird.c %>% filter(species == 'NOCA')
sosp.c <- bird.c %>% filter(species == 'SOSP')


# Process the data

amro.t.proc <- process.data(amro.t, groups = "sex")
cach.t.proc <- process.data(cach.t, groups = "sex")
carw.t.proc <- process.data(carw.t, groups = "sex")
grca.t.proc <- process.data(grca.t, groups = "sex")
howr.t.proc <- process.data(howr.t, groups = "sex")
noca.t.proc <- process.data(noca.t, groups = "sex")
sosp.t.proc <- process.data(sosp.t, groups = "sex")

amro.c.proc <- process.data(amro.c, groups = "sex")
cach.c.proc <- process.data(cach.c, groups = "sex")
carw.c.proc <- process.data(carw.c, groups = "sex")
grca.c.proc <- process.data(grca.c, groups = "sex")
howr.c.proc <- process.data(howr.c, groups = "sex")
noca.c.proc <- process.data(noca.c, groups = "sex")
sosp.c.proc <- process.data(sosp.c, groups = "sex")


# Design parameters
design.Phi <- list(static=c('imp','sex','bci','imp2', 'can', 'can2', 'tCats','cCats'))
design.p <- list(static = c('sex'))
design.parameters <- list(Phi = design.Phi, p = design.p)

amro.t.ddl <- make.design.data(amro.t.proc, parameters = design.parameters)
cach.t.ddl <- make.design.data(cach.t.proc, parameters = design.parameters)
carw.t.ddl <- make.design.data(carw.t.proc, parameters = design.parameters)
grca.t.ddl <- make.design.data(grca.t.proc, parameters = design.parameters)
howr.t.ddl <- make.design.data(howr.t.proc, parameters = design.parameters)
noca.t.ddl <- make.design.data(noca.t.proc, parameters = design.parameters)
sosp.t.ddl <- make.design.data(sosp.t.proc, parameters = design.parameters)

amro.c.ddl <- make.design.data(amro.c.proc, parameters = design.parameters)
cach.c.ddl <- make.design.data(cach.c.proc, parameters = design.parameters)
carw.c.ddl <- make.design.data(carw.c.proc, parameters = design.parameters)
grca.c.ddl <- make.design.data(grca.c.proc, parameters = design.parameters)
howr.c.ddl <- make.design.data(howr.c.proc, parameters = design.parameters)
noca.c.ddl <- make.design.data(noca.c.proc, parameters = design.parameters)
sosp.c.ddl <- make.design.data(sosp.c.proc, parameters = design.parameters)



# Functions that fit fourteen different models to transect, camera data

fit.models.t <- function(proc,ddl){
  Phi.dot <- list(formula = ~1)
  
  Phi.one <- list(formula = ~sex)
  Phi.two <- list(formula = ~sex*bci)
  Phi.three <- list(formula = ~sex*bci + imp)
  Phi.four <- list(formula = ~sex*bci + imp + imp2)
  Phi.five <- list(formula = ~sex*bci + imp*sex)
  Phi.six <- list(formula = ~sex*bci + imp*sex + imp2)
  
  Phi.seven <- list(formula= ~sex*tCats)
  
  Phi.eight <- list(
    formula= ~sex*bci + sex*tCats)
  
  Phi.nine <- list(
    formula= ~sex*bci + sex*tCats + bci*tCats)
  
  Phi.ten <- list(
    formula= ~sex*bci + sex*tCats + bci*tCats + imp*tCats)
  
  Phi.eleven <- list(
    formula= ~sex*bci + sex*tCats + bci*tCats + imp*tCats + imp2)
  
  Phi.twelve <- list(
    formula= ~sex*bci + sex*tCats + bci*tCats + imp*tCats + imp*sex)
  
  Phi.thirteen <- list(
    formula= ~sex*bci + sex*tCats + bci*tCats + imp*tCats + imp*sex + imp2)
  
  p.sex <- list(formula = ~sex)
  
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml,data=proc,ddl=ddl,external=FALSE,accumulate=FALSE)
  return(results)
}


fit.models.c <- function(proc,ddl){
  Phi.dot <- list(formula = ~1)
  
  Phi.one <- list(formula = ~sex)
  Phi.two <- list(formula = ~sex*bci)
  Phi.three <- list(formula = ~sex*bci + imp)
  Phi.four <- list(formula = ~sex*bci + imp + imp2)
  Phi.five <- list(formula = ~sex*bci + imp*sex)
  Phi.six <- list(formula = ~sex*bci + imp*sex + imp2)
  
  Phi.seven <- list(formula= ~sex*cCats)
  
  Phi.eight <- list(
    formula= ~sex*bci + sex*cCats)
  
  Phi.nine <- list(
    formula= ~sex*bci + sex*cCats + bci*cCats)
  
  Phi.ten <- list(
    formula= ~sex*bci + sex*cCats + bci*cCats + imp*cCats)
  
  Phi.eleven <- list(
    formula= ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp2)
  
  Phi.twelve <- list(
    formula= ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp*sex)
  
  Phi.thirteen <- list(
    formula= ~sex*bci + sex*cCats + bci*cCats + imp*cCats + imp*sex + imp2)
  
  p.sex <- list(formula = ~sex)
  
  cml <- create.model.list(c("Phi","p"))
  results <- crm.wrapper(cml,data=proc,ddl=ddl,external=FALSE,accumulate=FALSE)
  return(results)
}


# Fit the models

amro.models.t <- fit.models.t(amro.t.proc,amro.t.ddl)
cach.models.t <- fit.models.t(cach.t.proc,cach.t.ddl)
carw.models.t <- fit.models.t(carw.t.proc,carw.t.ddl)
grca.models.t <- fit.models.t(grca.t.proc,grca.t.ddl)
howr.models.t <- fit.models.t(howr.t.proc,howr.t.ddl)
noca.models.t <- fit.models.t(noca.t.proc,noca.t.ddl)
sosp.models.t <- fit.models.t(sosp.t.proc,sosp.t.ddl)

amro.models.c <- fit.models.c(amro.c.proc,amro.c.ddl)
cach.models.c <- fit.models.c(cach.c.proc,cach.c.ddl)
carw.models.c <- fit.models.c(carw.c.proc,carw.c.ddl)
grca.models.c <- fit.models.c(grca.c.proc,grca.c.ddl)
howr.models.c <- fit.models.c(howr.c.proc,howr.c.ddl)
noca.models.c <- fit.models.c(noca.c.proc,noca.c.ddl)
sosp.models.c <- fit.models.c(sosp.c.proc,sosp.c.ddl)


