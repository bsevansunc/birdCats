catSites <- read.csv('catSites.csv',
                     stringsAsFactors = FALSE) %>%
  tbl_df %>%
  .$site

temporalSamplingFun <- function(){
  dfList <- vector('list', length = 3)
  for(i in 1:3){
    dfList[[i]] <- data.frame(
      site = catSites,
      monthSample = rep(i, 52),
      eventSample = sample(1:52, 52)
    ) 
  }
  df <- do.call('rbind', dfList) %>% 
    tbl_df %>% 
    arrange(monthSample, eventSample) %>%
    mutate(eventWk =  reshape::round_any(eventSample/13, 1, ceiling)) %>%
    select(site, monthSample, eventWk, eventSample)
  return(df)
}


test1 <- df %>%
  group_by(site) %>%
  mutate(t1 = eventWk == 4 & lead(eventWk == 1)) %>%
  filter(t1 == TRUE)




test2 <- temporalSamplingFun() %>%
  group_by(site) %>%
  mutate(t1 = eventWk == 4 & lead(eventWk == 1)) %>%
  filter(t1 == TRUE)

nrowVector = numeric()

for(i in 1:10000){
  nrowVector[i] <- nrow(
    temporalSamplingFun() %>%
    group_by(site) %>%
    mutate(t1 = eventWk == 4 & lead(eventWk == 1)) %>%
    filter(t1 == TRUE)
  )
}

dfList = list()

for(i in 1:10000){
  df <- temporalSamplingFun()
  if(nrow(
    df %>%
    group_by(site) %>%
    mutate(t1 = eventWk == 4 & lead(eventWk == 1)) %>%
    filter(t1 == TRUE)
  ) < 1) {
    dfList[[i]] <- df %>%
      mutate(iteration = i)
  } else {
    dfList[[i]] <- NULL
    }
}

randomizedSamples <- bind_rows(dfList)


