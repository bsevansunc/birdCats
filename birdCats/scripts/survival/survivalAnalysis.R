library('tidyverse')
library('marked')
library('stringr')

# Read in data

data <- read_csv('data/ch.csv')


# Filter out individuals of unknown sex and body condition, add imp^2
bird <- data %>%
  mutate(imp2 = imp^2) %>%
  filter(!is.na(bci)) %>%
  filter(sex != 'U') %>%
  as.data.frame


# Set sex as a factor with levels 'M' and 'F'
bird$sex <- as.factor(bird$sex)


# By species
amro <- bird %>% filter(speciesEnc == 'AMRO')
cach <- bird %>% filter(speciesEnc == 'CACH')
carw <- bird %>% filter(speciesEnc == 'CARW')
grca <- bird %>% filter(speciesEnc == 'GRCA')
howr <- bird %>% filter(speciesEnc == 'HOWR')
noca <- bird %>% filter(speciesEnc == 'NOCA')
nomo <- bird %>% filter(speciesEnc == 'NOMO')
sosp <- bird %>% filter(speciesEnc == 'SOSP')


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
Phi.noCats <- list(formula = ~imp+sex+bci+sex*bci+imp2)
Phi.tCats <- list(formula = ~imp+tCats+sex+bci+bci*tCats+sex*tCats+sex*bci+imp2)
Phi.cCats <- list(formula = ~imp+cCats+sex+bci+bci*cCats+sex*cCats+sex*bci+imp2)
Phi.mCats <- list(formula = ~imp+mCats+sex+bci+bci*mCats+sex*mCats+sex*bci+imp2)

p.dot <- list(formula = ~1)
p.sex <- list(formula = ~sex)
p.sa <- list(formula = ~sex+active)

# Models
# Combined species
all.null <- crm(data=all.proc, ddl=all.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sa))

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
  Phi = Phi.dot, p = p.sex))

noca.cCats <- crm(data=noca.proc, ddl=noca.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# CACH
cach.null <- crm(data=cach.proc, ddl=cach.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

cach.cCats <- crm(data=cach.proc, ddl=cach.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# CARW
carw.null <- crm(data=carw.proc, ddl=carw.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

carw.cCats <- crm(data=carw.proc, ddl=carw.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# GRCA
grca.null <- crm(data=grca.proc, ddl=grca.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

grca.cCats <- crm(data=grca.proc, ddl=grca.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# HOWR
howr.null <- crm(data=howr.proc, ddl=howr.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

howr.cCats <- crm(data=howr.proc, ddl=howr.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# SOSP
sosp.null <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

sosp.cCats <- crm(data=sosp.proc, ddl=sosp.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# AMRO
amro.null <- crm(data=amro.proc, ddl=amro.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

amro.cCats <- crm(data=amro.proc, ddl=amro.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))

# NOMO
nomo.null <- crm(data=nomo.proc, ddl=nomo.ddl, model.parameters = list(
  Phi = Phi.dot, p = p.sex))

nomo.cCats <- crm(data=nomo.proc, ddl=nomo.ddl, model.parameters=list(
  Phi=Phi.cCats, p=p.sex))



