cAbund <- cAbund[,-1]
tAbund <- tAbund[,-1]

covs <- cAbund[,7:13]


arrays <- vector(mode = 'numeric', length = 10)

carrays <- vector(mode = 'numeric', length = 10)


cor.array <- function(data){
  tarray <- vector(mode = 'numeric', length = 10)
  for(i in 1:length(names(data))){
    tarray[i]
  }
}

arrays[1] <- cor.test(tAbund$hDensity,tAbund$medianIncome, method = 'spearman')$p.value
arrays[2] <- cor.test(tAbund$hDensity,tAbund$age, method = 'spearman')$p.value
arrays[3] <- cor.test(tAbund$hDensity,tAbund$marred, method = 'spearman')$p.value
arrays[4] <- cor.test(tAbund$hDensity,tAbund$eduHS, method = 'spearman')$p.value

arrays[5] <- cor.test(tAbund$medianIncome,age, method = 'spearman')$p.value
arrays[6] <- cor.test(tAbund$medianIncome,tAbund$marred, method = 'spearman')$p.value
arrays[7] <- cor.test(tAbund$medianIncome,tAbund$eduHS, method = 'spearman')$p.value

arrays[8] <- cor.test(tAbund$age,tAbund$marred, method = 'spearman')$p.value
arrays[9] <- cor.test(tAbund$age,tAbund$eduHS, method = 'spearman')$p.value

arrays[10] <- cor.test(tAbund$marred,eduHS, method = 'spearman')$p.value



carrays[1] <- cor.test(cAbund$hDensity,cAbund$medianIncome, method = 'spearman')$p.value
carrays[2] <- cor.test(cAbund$hDensity,cAbund$age, method = 'spearman')$p.value
carrays[3] <- cor.test(cAbund$hDensity,cAbund$marred, method = 'spearman')$p.value
carrays[4] <- cor.test(cAbund$hDensity,cAbund$eduHS, method = 'spearman')$p.value

carrays[5] <- cor.test(cAbund$medianIncome,age, method = 'spearman')$p.value
carrays[6] <- cor.test(cAbund$medianIncome,cAbund$marred, method = 'spearman')$p.value
carrays[7] <- cor.test(cAbund$medianIncome,cAbund$eduHS, method = 'spearman')$p.value

carrays[8] <- cor.test(cAbund$age,cAbund$marred, method = 'spearman')$p.value
carrays[9] <- cor.test(cAbund$age,cAbund$eduHS, method = 'spearman')$p.value

carrays[10] <- cor.test(cAbund$marred,cAbund$eduHS, method = 'spearman')$p.value