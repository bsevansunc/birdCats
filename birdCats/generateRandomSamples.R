# Read sites for study:

catSites <- read.csv('catSites.csv',
                     stringsAsFactors = FALSE) %>%
  tbl_df %>%
  .$site

# Function generates a random sampling schedule for sites, returns data frame

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


dfList = vector('list', length = 10000)

# For loop checks whether sites are visited on the 4th then 1st week. If this is the case, it returns NULL, if this is not the case, it returns the randomly sampled data frame:

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

# Binds elements in the rows together into a single data frame then selects the first iteration that matches the 4 to 1 condition.

randomizedSamples <- bind_rows(dfList) %>%
  filter(iteration == min(iteration))


