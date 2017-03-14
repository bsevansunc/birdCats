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
design.Phi <- list(static = c('imp','sex','bci','tCats','cCats','mCats','imp2'))
design.p <- list(static = c('sex'), time.varying = c("active"))
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
Phi.dot <- list(formula = ~1)
Phi.cCats <-list(formula = ~imp+cCats)

Phi.one <- list(formula = ~sex)
Phi.two <- list(formula = ~bci)
Phi.three <- list(formula = ~sex+bci)
Phi.four <- list(formula = ~sex+bci+sex*bci)
Phi.five <- list(formula = ~sex+bci+imp)
Phi.six <- list(formula = ~sex+bci+imp+imp2)
Phi.seven <- list(formula = ~sex+bci+sex*bci+imp)
Phi.eight <- list(formula = ~sex+bci+sex*bci+imp+imp2)

Phi.nine <- list(formula = ~sex+bci+cCats)
Phi.ten <- list(formula = ~sex+bci+cCats+cCats*sex)
Phi.eleven <- list(formula = ~sex+bci+cCats+cCats*bci)
Phi.twelve <- list(formula = ~sex+bci+cCats+sex*cCats+bci*cCats)
Phi.thirteen <- list(formula = ~imp+sex+bci+cCats+sex*cCats+bci*cCats)

Phi.tCats <- list(formula = ~imp+tCats+sex+bci+bci*tCats+sex*tCats+sex*bci+imp2)
Phi.cCats <- list(formula = ~imp+cCats+sex+bci+bci*cCats+sex*cCats+sex*bci+imp2)
Phi.mCats <- list(formula = ~imp+mCats+sex+bci+bci*mCats+sex*mCats+sex*bci+imp2)

p.dot <- list(formula = ~1)
p.sex <- list(formula = ~sex)
p.active <- list(formula = ~active)
p.sa <- list(formula = ~sex+active)


# Models
# Combined species
all.null <- crm(data=all.proc, ddl=all.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))


all.one <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.one, p=p.sa))

all.two <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.two, p=p.sa))

all.three <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.three, p=p.sa))

all.four <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.four, p=p.sa))

all.five <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.five, p=p.sa))

all.six <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.six, p=p.sa))

all.seven <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.seven, p=p.sa))

all.eight <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.eight, p=p.sa))

all.nine <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))

all.ten <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.ten, p=p.sa))

all.eleven <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.eleven, p=p.sa))

all.twelve <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.twelve, p=p.sa))

all.thirteen <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.thirteen, p=p.sa))



all.tCats <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.tCats, p=p.sa))

all.cCats <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

all.mCats <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.mCats, p=p.sa))

all.noCats <- crm(data=all.proc, ddl=all.ddl, model.parameters=list(
  Phi=Phi.noCats, p=p.sa))





# NOCA
noca.null <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

noca.one <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.one, p = p.sa))

noca.two <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.two, p = p.sa))

noca.three <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.three, p = p.sa))

noca.nine <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.nine, p = p.sa))

noca.twelve <- crm(data=noca.proc, ddl=noca.ddl, model.parameters = list(
  Phi = Phi.twelve, p = p.sa))

noca.cCats <- crm(data=noca.proc, ddl=noca.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

# CACH
cach.null <- crm(data=cach.proc, ddl=cach.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

cach.cCats <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

cach.one <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.one, p=p.sa))

cach.two <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.two, p=p.sa))

cach.three <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.three, p=p.sa))

cach.four <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.four, p=p.sa))

cach.five <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.five, p=p.sa))

cach.six <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.six, p=p.sa))

cach.seven <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.seven, p=p.sa))

cach.eight <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.eight, p=p.sa))

cach.nine <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))

cach.ten <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.ten, p=p.sa))

cach.eleven <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.eleven, p=p.sa))

cach.twelve <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.twelve, p=p.sa))

cach.thirteen <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.thirteen, p=p.sa))

# CARW
carw.null <- crm(data=carw.proc, ddl=carw.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

carw.cCats <- crm(data=carw.proc, ddl=carw.ddl, model.parameters = list(
  Phi = Phi.cCats, p = p.sa))

carw.one <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.one, p=p.sa))

carw.two <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.two, p=p.sa))

carw.three <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.three, p=p.sa))

carw.four <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.four, p=p.sa))

carw.five <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.five, p=p.sa))

carw.six <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.six, p=p.sa))

carw.seven <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.seven, p=p.sa))

carw.eight <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.eight, p=p.sa))

carw.nine <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))

carw.ten <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.ten, p=p.sa))

carw.eleven <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.eleven, p=p.sa))

carw.twelve <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.twelve, p=p.sa))

carw.thirteen <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.thirteen, p=p.sa))

# GRCA
grca.null <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

grca.one <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.one, p = p.sa))

grca.two <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.two, p = p.sa))

grca.three <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.three, p = p.sa))

grca.eight <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.eight, p = p.sa))

grca.nine <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.nine, p = p.sa))

grca.ten <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.ten, p = p.sa))

grca.eleven <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.eleven, p = p.sa))

grca.twelve <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.twelve, p = p.sa))

grca.thirteen <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.thirteen, p = p.sa))

grca.cCats <- crm(data=grca.proc, ddl=grca.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

# HOWR
howr.null <- crm(data=howr.proc, ddl=howr.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

howr.cCats <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

howr.three <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.three, p=p.sa))

howr.eight <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.eight, p=p.sa))

howr.nine <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))

howr.twelve <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.twelve, p=p.sa))

# SOSP
sosp.null <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

sosp.one <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.one, p = p.sa))

sosp.two <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.two, p = p.sa))

sosp.three <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.three, p = p.sa))

sosp.four <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.four, p = p.sa))

sosp.eight <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.eight, p = p.sa))

sosp.nine <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.nine, p = p.sa))

sosp.ten <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.ten, p = p.sa))

sosp.eleven <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.eleven, p = p.sa))

sosp.twelve <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.twelve, p = p.sa))

sosp.thirteen <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.thirteen, p = p.sa))

sosp.cCats <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

# AMRO
amro.null <- crm(data=amro.proc, ddl=amro.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

amro.cCats <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

amro.one <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.one, p=p.sa))

amro.two <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.two, p=p.sa))

amro.three <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.three, p=p.sa))

amro.nine <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))

# NOMO
nomo.null <- crm(data=nomo.proc, ddl=nomo.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

nomo.cCats <- crm(data=nomo.proc, ddl=nomo.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sa))

nomo.nine <- crm(data=nomo.proc, ddl=nomo.ddl, model.parameters=list(
  Phi=Phi.nine, p=p.sa))



