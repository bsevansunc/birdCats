library('tidyverse')
library('marked')
library('stringr')

removeSites <- c('GERYERIMD1', 'OLONMARDC1', 'MISSEDDC1', 'WOLFKARDC1', 'WOLFAMYDC1')

# Read in data

data <- read_csv('data/ch.csv')


# Filter out individuals of unknown sex and body condition, add imp^2
bird <- data %>%
  mutate(imp2 = imp^2) %>%
  filter(
    !is.na(bci),
    age != "U" & age != 'HY',
    sex != 'U'
#    ,!site %in% removeSites
    ) %>%
  as.data.frame


# Set sex as a factor with levels 'M' and 'F'
bird$sex <- as.factor(bird$sex)


# By species
amro <- bird %>% filter(species == 'AMRO')
cach <- bird %>% filter(species == 'CACH')
carw <- bird %>% filter(species == 'CARW')
grca <- bird %>% filter(species == 'GRCA')
howr <- bird %>% filter(species == 'HOWR')
noca <- bird %>% filter(species == 'NOCA')
nomo <- bird %>% filter(species == 'NOMO')
sosp <- bird %>% filter(species == 'SOSP')


# Process the data
all.proc <- process.data(bird, groups = "sex")

amro.proc <- process.data(amro, groups = "sex")
cach.proc <- process.data(cach, groups = "sex")
carw.proc <- process.data(carw, groups = "sex")
grca.proc <- process.data(grca, groups = "sex")
howr.proc <- process.data(howr, groups = "sex")
noca.proc <- process.data(noca, groups = "sex")
nomo.proc <- process.data(nomo, groups = "sex")
sosp.proc <- process.data(sosp, groups = "sex")


# Design parameters
design.Phi <- list(static = c('imp','sex','bci','imp2', 'can', 'can2', 'avgTrans'))
design.p <- list(static = c('sex'))
design.parameters <- list(Phi = design.Phi, p = design.p)

all.ddl <- make.design.data(all.proc, parameters = design.parameters)

amro.ddl <- make.design.data(amro.proc, parameters = design.parameters)
cach.ddl <- make.design.data(cach.proc, parameters = design.parameters)
carw.ddl <- make.design.data(carw.proc, parameters = design.parameters)
grca.ddl <- make.design.data(grca.proc, parameters = design.parameters)
howr.ddl <- make.design.data(howr.proc, parameters = design.parameters)
noca.ddl <- make.design.data(noca.proc, parameters = design.parameters)
nomo.ddl <- make.design.data(nomo.proc, parameters = design.parameters)
sosp.ddl <- make.design.data(sosp.proc, parameters = design.parameters)




# Model parameters

fit.models <- function(proc,ddl){
  
Phi.dot <- list(formula = ~1)
Phi.trans <- list(formula = ~avgTrans)
Phi.sexTrans <- list(formula = ~sex+avgTrans)

Phi.one <- list(formula = ~sex)
Phi.two <- list(formula = ~bci)
Phi.three <- list(formula = ~sex+bci)
Phi.four <- list(formula = ~sex+bci+sex*bci)
Phi.five <- list(formula = ~sex+bci+imp)
Phi.five.b <- list(formula = ~sex+bci+can)
Phi.six <- list(formula = ~sex+bci+imp+imp2)
Phi.six.b <- list(formula = ~sex+bci+can+can2)
Phi.seven <- list(formula = ~sex+bci+sex*bci+imp)
Phi.seven.b <- list(formula = ~sex+bci+sex*bci+can)
Phi.eight <- list(formula = ~sex+bci+sex*bci+imp+imp2)
Phi.eight.b <- list(formula = ~sex+bci+sex*bci+can+can2)

Phi.nine <- list(formula = ~sex+bci+avgTrans)
Phi.ten <- list(formula = ~sex+bci+avgTrans+sex*avgTrans)
Phi.eleven <- list(formula = ~sex+bci+avgTrans+bci*avgTrans)
Phi.twelve <- list(formula = ~sex+bci+avgTrans+sex*avgTrans+bci*avgTrans)
Phi.thirteen <- list(formula = ~imp+sex+bci+avgTrans+sex*avgTrans+bci*avgTrans)
Phi.thirteen.b <- list(formula = ~can+sex+bci+avgTrans+sex*avgTrans+bci*avgTrans)
Phi.fourteen <- list(formula = ~imp+sex+bci+avgTrans+sex*avgTrans+bci*avgTrans+imp2)
Phi.fourteen.b <- list(formula = ~can+sex+bci+avgTrans+sex*avgTrans+bci*avgTrans+can2)


p.sex <- list(formula = ~sex)

cml <- create.model.list(c("Phi","p"))

results <- crm.wrapper(cml,data=proc,ddl=ddl,external=FALSE,accumulate=FALSE)

return(results)
}

amro.models <- fit.models(amro.proc,amro.ddl)
cach.models <- fit.models(cach.proc,cach.ddl)
carw.models <- fit.models(carw.proc,carw.ddl)
grca.models <- fit.models(grca.proc,grca.ddl)
howr.models <- fit.models(howr.proc,howr.ddl)
noca.models <- fit.models(noca.proc,noca.ddl)
sosp.models <- fit.models(sosp.proc,sosp.ddl)
nomo.models <- fit.models(nomo.proc,nomo.ddl)
