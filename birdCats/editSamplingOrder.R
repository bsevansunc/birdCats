# Working directory: Kevin's laptop: C:/Users/Kevin/Documents/gitHub/birdCats/birdCats

library(RCurl)
library(dplyr)
library(tidyr)


# Get sampling order:

samplingOrder <- read.csv('samplingOrder.csv', stringsAsFactors = FALSE) %>%
  tbl_df()

# ------------------------------------------------------------------------------------------*
# ---- CHANGE THE SAMPLING ORDER ----
# ==========================================================================================*


# Give GERYERIMD1 order to LEONTESDC1:

samplingOrder[samplingOrder$site == 'GERYERIMD1', 'site'] <-'LEONTESDC1a'  


# Give WOLFAMYDC1 order to PIROELADC1:

samplingOrder[samplingOrder$site == 'WOLFAMYDC1', 'site'] <-'PIROELADC1a'


# Give FREYFREMD1 order to BRIGKIMMD1:

samplingOrder[samplingOrder$site == 'FREYFREMD1', 'site'] <-'BRIGKIMMD1a'


# Give LEONTESDC1 order to FISHFERVA1:

samplingOrder[samplingOrder$site == 'LEONTESDC1', 'site'] <-'FISHFERVA1a'


# Give FISHFERVA1 order to GERYERIMD1:

samplingOrder[samplingOrder$site == 'FISHFERVA1', 'site'] <-'GERYERIMD1a'


# Give PIROELADC1 order to WOLFAMYDC1:

samplingOrder[samplingOrder$site == 'PIROELADC1', 'site'] <-'WOLFAMYDC1a'


# Give BRIGKIMMD1 order to FREYFREMD1:

samplingOrder[samplingOrder$site == 'BRIGKIMMD1', 'site'] <-'FREYFREMD1a'


# ---------------------------------------------------------------------------------------*
# ---- DEAL WITH THE MESS: ----
# =======================================================================================*



samplingOrder[samplingOrder$site == 'GERYERIMD1a', 'site'] <-'GERYERIMD1'

samplingOrder[samplingOrder$site == 'LEONTESDC1a', 'site'] <-'LEONTESDC1'

samplingOrder[samplingOrder$site == 'FISHFERVA1a', 'site'] <-'FISHFERVA1'

samplingOrder[samplingOrder$site == 'WOLFAMYDC1a', 'site'] <-'WOLFAMYDC1'

samplingOrder[samplingOrder$site == 'PIROELADC1a', 'site'] <-'PIROELADC1'

samplingOrder[samplingOrder$site == 'FREYFREMD1a', 'site'] <-'FREYFREMD1'

samplingOrder[samplingOrder$site == 'BRIGKIMMD1a', 'site'] <-'BRIGKIMMD1'


# ---------------------------------------------------------------------------------------*
# ---- SAVE THE FILE ----
#========================================================================================*



write.csv(samplingOrder, file = 'samplingOrder.csv')
