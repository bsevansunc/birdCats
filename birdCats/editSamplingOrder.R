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


# Give LEONTESDC1 order to FISHFERVA1:

samplingOrder[samplingOrder$site == 'LEONTESDC1', 'site'] <-'FISHFERVA1a'


# Give FISHFERVA1 order to GERYERIMD1:

samplingOrder[samplingOrder$site == 'FISHFERVA1', 'site'] <-'GERYERIMD1a'


# Give WOLFAMYDC1 order to PIROELADC1 and vice versa:

samplingOrder[samplingOrder$site == 'WOLFAMYDC1', 'site'] <-'PIROELADC1a'
samplingOrder[samplingOrder$site == 'PIROELADC1', 'site'] <-'WOLFAMYDC1a'


# Give FREYFREMD1 order to BRIGKIMMD1 and vice versa:

samplingOrder[samplingOrder$site == 'FREYFREMD1', 'site'] <-'BRIGKIMMD1a'
samplingOrder[samplingOrder$site == 'BRIGKIMMD1', 'site'] <-'FREYFREMD1a'


# Give BROWANNMD1 order to BROWANNDC1 and vice versa:

samplingOrder[samplingOrder$site == 'BROWANNMD1', 'site'] <-'BROWANNDC1a'
samplingOrder[samplingOrder$site == 'BROWANNDC1', 'site'] <-'BROWANNMD1a'


# Replace BOONGREMD1 with STUECHRVA2:

samplingOrder[samplingOrder$site == 'BOONGREMD1', 'site'] <-'STUECHRVA2'


# Replace PERNMIKMD1 with ROHRSALMD1:

samplingOrder[samplingOrder$site == 'PERNMIKMD1', 'site'] <-'ROHRSALMD1'


# Replace FREYFREMD1 with STANJEFMD1:

samplingOrder[samplingOrder$site == 'FREYFREMD1', 'site'] <-'STANJEFMD1'



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

samplingOrder[samplingOrder$site == 'BROWANNDC1a', 'site'] <-'BROWANNDC1'

samplingOrder[samplingOrder$site == 'BROWANNMD1a', 'site'] <-'BROWANNMD1'



# ---------------------------------------------------------------------------------------*
# ---- SAVE THE FILE ----
#========================================================================================*



write.csv(samplingOrder, file = 'samplingOrder.csv', row.names = FALSE)
