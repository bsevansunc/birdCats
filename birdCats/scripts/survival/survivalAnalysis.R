#--------------------------------------------------------------------*
# ---- Set up ----
#--------------------------------------------------------------------*

library('tidyverse')
library('marked')
library('stringr')
library('AICcmodavg')
library('MuMIn')

removeSites <- c('GERYERIMD1', 'OLONMARDC1', 'MISSEDDC1',
                 'WOLFKARDC1', 'WOLFAMYDC1')

# Read and format data:

data <- read_csv('data/ch.csv') %>%
  filter(
    age %in% c('AHY', 'SY', 'ASY'),
    sex %in% c('M','F'),
    !is.na(bci)) %>%
  as.data.frame() %>%
  mutate(sex = as.factor(sex),
         cCats = cCats*3,
         tCats = tCats*6) %>%
  select(ch, bandNumber,site, bandYear, age, sex, species, bci:tCats)



#--------------------------------------------------------------------*
# ---- Functions ----
#--------------------------------------------------------------------*

# Function to scale covariates:

scaleCovs <- function(data){
  imp <- data %>% select(site,imp) %>% distinct()
  imp$imp <- scale(imp$imp)[,1]
  data <- data %>% select(-imp)
  data <- left_join(data, imp, by = 'site')
  
  imp2 <- data %>% select(site,imp2) %>% distinct()
  imp2$imp2 <- scale(imp2$imp2)[,1]
  data <- data %>% select(-imp2)
  data <- left_join(data, imp2, by = 'site')
  
  can <- data %>% select(site,can) %>% distinct()
  can$can <- scale(can$can)[,1]
  data <- data %>% select(-can)
  data <- left_join(data, can, by = 'site')
  
  can2 <- data %>% select(site,can2) %>% distinct()
  can2$can2 <- scale(can2$can2)[,1]
  data <- data %>% select(-can2)
  data <- left_join(data, can2, by = 'site')
  
  cats <- data %>% select(site,cat) %>% distinct()
  cats$cat <- scale(cats$cat)[,1]
  data <- data %>% select(-cat)
  data <- left_join(data, cats, by = 'site')
  
  return(data)
}

# Function to get frame for a given species:
# Issue with this that I'll fix later--if you scale inside the full
# frame, I think the means and stdevs will be off.

getSppFrame <- function(data, spp, catVar, scaled = TRUE){
  sppData <- data %>%
    filter(species == spp) %>%
    mutate_(cat = catVar) %>%
    select(-c(cCats,tCats))
  if(scaled == TRUE){
    sppData <- scaleCovs(sppData)
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
                         external=FALSE,accumulate=FALSE,hessian=TRUE)
  return(results)
}


# Function that gives an AICc table, rather than an AIC table

AICc.tab <- function(modelList,spp){
  LL <- model.table(modelList[[spp]])[,6]*-0.5
  ss <- data %>% filter(species == spp) %>% nrow
  K <- model.table(modelList[[spp]])[,2]
  names <- model.table(modelList[[spp]])[,1]
  
  aictab <- aictabCustom(LL,modnames=names,K,nobs=ss)
  return(aictab)
}


# Function that calculates the weighted betas for a given sp and cov

weightedBeta <- function(modellist, tablelist, spp, cov){
  betas <- vector('numeric',length=76)
  for(i in 1:76){
    betas[i] <- modellist[[spp]][[i]]$results$beta$Phi[cov]
  }
  LLs <- vector('numeric',length=76)
  for(i in 1:76){
    LLs[i] <- (modellist[[spp]][[i]]$results$neg2lnl)*-0.5
  }
  
  SEs <- vector('numeric',length=76)
  for(i in 1:76){
    SEs[i] <- coef(modellist[[spp]][[i]])[paste('Phi.',cov,sep=''),'se']
  }
  
  betadf <- data.frame(beta=betas,LL=LLs,SE=SEs)
  betadf <- betadf[order(betadf$LL),]
  
  table <- tablelist[[spp]]
  table <- table[order(table$LL),]
  
  wtB <- sum(na.omit(table$AICcWt*betadf$beta)[1:length(na.omit(betadf$beta))])
  wtSE <- sum(na.omit(table$AICcWt*betadf$SE)[1:length(na.omit(betadf$SE))])
  
  
  return(paste(wtB,'+/-',wtSE, sep = ' '))
}


#--------------------------------------------------------------------*
# ---- Model running ----
#--------------------------------------------------------------------*

# Fit the models:

#data <- scaleCovs(data)

sppVector <- data$species %>% unique %>% sort

modelListSppCamera <- vector('list', length = length(sppVector))
modelListSppTransect <- vector('list', length = length(sppVector))

modelList <- vector('list', length = length(sppVector))

names(modelListSppCamera) <- sppVector
names(modelListSppTransect) <- sppVector
names(modelList) <- sppVector



for(i in 1:length(sppVector)){
  modelListSppCamera[[i]] <- fitModels(
    data = data %>% filter(!site %in% removeSites),
    spp = sppVector[i],
    catVar = 'cCats',
    scaled = TRUE)
  modelListSppTransect[[i]] <- fitModels(
    data = data,
    spp = sppVector[i],
    catVar = 'tCats',
    scaled = TRUE)
}



# Example output:

modelList$GRCA

modelList$GRCA$Phi.21.p.dot

modelList$SOSP$Phi.7.p.sex

AMROtab <- AICc.tab(modelList,'AMRO')
CACHtab <- AICc.tab(modelList,'CACH')
CARWtab <- AICc.tab(modelList,'CARW')
GRCAtab <- AICc.tab(modelList,'GRCA')
HOWRtab <- AICc.tab(modelList,'HOWR')
NOCAtab <- AICc.tab(modelList,'NOCA')
SOSPtab <- AICc.tab(modelList,'SOSP')

tableList <- list(AMROtab,CACHtab,CARWtab,GRCAtab,HOWRtab,NOCAtab,SOSPtab)
names(tableList) <- c('AMRO','CACH','CARW','GRCA','HOWR','NOCA','SOSP')

AMROtabT <- AICc.tab(modelListSppTransect,'AMRO')
CACHtabT <- AICc.tab(modelListSppTransect,'CACH')
CARWtabT <- AICc.tab(modelListSppTransect,'CARW')
GRCAtabT <- AICc.tab(modelListSppTransect,'GRCA')
HOWRtabT <- AICc.tab(modelListSppTransect,'HOWR')
NOCAtabT <- AICc.tab(modelListSppTransect,'NOCA')
SOSPtabT <- AICc.tab(modelListSppTransect,'SOSP')

tableListT <- list(AMROtabT,CACHtabT,CARWtabT,GRCAtabT,HOWRtabT,NOCAtabT,SOSPtabT)
names(tableListT) <- c('AMRO','CACH','CARW','GRCA','HOWR','NOCA','SOSP')

AMROtabC <- AICc.tab(modelListSppCamera,'AMRO')
CACHtabC <- AICc.tab(modelListSppCamera,'CACH')
CARWtabC <- AICc.tab(modelListSppCamera,'CARW')
GRCAtabC <- AICc.tab(modelListSppCamera,'GRCA')
HOWRtabC <- AICc.tab(modelListSppCamera,'HOWR')
NOCAtabC <- AICc.tab(modelListSppCamera,'NOCA')
SOSPtabC <- AICc.tab(modelListSppCamera,'SOSP')


tableListC <- list(AMROtab,CACHtab,CARWtab,GRCAtab,HOWRtab,NOCAtab,SOSPtab)
names(tableListC) <- c('AMRO','CACH','CARW','GRCA','HOWR','NOCA','SOSP')


modelList$SOSP$Phi.7.p.sex$results$reals$Phi %>%
  ggplot(aes(x = cat, y = estimate)) +
  geom_point(cex = 3, alpha = .5) +
  theme_bw()

modelList$GRCA$Phi.7.p.sex$results$reals$Phi %>%
  ggplot(aes(x = cat, y = estimate)) +
  geom_point(cex = 3, alpha = .5) +
  theme_bw()


weightedBeta(modelList, tableList, 'CARW', 'can')






