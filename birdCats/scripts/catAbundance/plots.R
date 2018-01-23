# ==================================================================================*
# ---------------------------------- Set up ----------------------------------------
# ==================================================================================*

# Function searches packages in installed package list,
# add them if you don't have them, and loads the library:

smartLibrary <- function(packageVector){
  for(i in 1:length(packageVector)){
    package <- packageVector[i]
    if(!package %in% rownames(installed.packages())){
      install.packages(packageVector[i],repos="http://cran.rstudio.com/",
                       dependencies=TRUE)
    }
  }
  lapply(packageVector, library, character.only = TRUE)
}

smartLibrary(
  c('dplyr', 'tidyr', 'ggplot2', 'AICcmodavg','MuMIn', 'lubridate')
)


options(stringsAsFactors=F)



# Read in the data

tAbund <- read.csv('data/abundanceTrans.csv')
cAbund <- read.csv('data/abundanceCam.csv')

NNSites <- read.csv('data/lc100.csv')

eduData <- read.csv('data/eduHS.csv')

covs <- read.csv('data/covariateData.csv') %>%
  select(-c(eduHS, eduC)) %>%
  left_join(eduData) %>%
  tbl_df %>%
  arrange(site)

camDet <- read.csv('data/catDataCamera.csv') %>%
  filter(species == 'cat') %>%
  mutate(time = hour(hm(time))*60 + minute(hm(time)))

transAll <- read.csv('data/catDataActivity.csv') %>%
  filter(activity == 'transect') %>%
  select(site, activity, time)

catSiteActivity <- read.csv('data/catDataActivity.csv') %>%
  tbl_df

transCat <- read.csv('data/catDataTransect.csv') %>%
  tbl_df %>%
  filter(!is.na(count)) %>%
  filter(species == 'cat', count != 0) %>%
  left_join(catSiteActivity %>%
              filter(activity == 'transect'),
            by = c('site', 'visit'))


# ==================================================================================*
# --------------------------- Plot abundance data ----------------------------------
# ==================================================================================*

# ----------------------------------------------------------------------------------*
# ---- Transect data plots ----
# ----------------------------------------------------------------------------------*


# Plot cat abundance by impervious surface (transect)

plotImpTrans <- ggplot(data=tAbund, mapping=aes(x=imp,y=abunds/19.774))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/19.774, ymax=(abunds+SE)/19.774))+
  scale_x_continuous("Impervious surface (%)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=10),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=10))


# Plot cat abundance by canopy cover (transect)

plotCanTrans <- ggplot(data=tAbund, mapping=aes(x=can,y=abunds/19.774))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/19.774, ymax=(abunds+SE)/19.774))+
  scale_x_continuous("Canopy cover (%)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=14),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=14))


# Plot cat abundance by human population density (transect)

plotDensTrans <- ggplot(data=tAbund, mapping=aes(x=hDensity,y=abunds/19.774))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/19.774, ymax=(abunds+SE)/19.774))+
  scale_x_continuous("Human population density (UNITS?)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=14),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=14))


# Plot cat abundance by high school education (transect)

plotEduTrans <- ggplot(data=tAbund, mapping=aes(x=eduHS*100,y=abunds/19.774))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/19.774, ymax=(abunds+SE)/19.774))+
  scale_x_continuous("High school degree (%)", 
                     breaks=c(80,85,90,95,100),
                     labels=c(80,85,90,95,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=10),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=10))


# ----------------------------------------------------------------------------------*
# ---- Camera data plots ----
# ----------------------------------------------------------------------------------*


# Plot cat abundance by impervious surface (camera)

plotImpCam <- ggplot(data=cAbund, mapping=aes(x=imp,y=abunds/1.197))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/1.197, ymax=(abunds+SE)/1.197))+
  scale_x_continuous("Impervious surface (%)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=10),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=10))



# Plot cat abundance by canopy cover (camera)

plotCanCam <- ggplot(data=cAbund, mapping=aes(x=can,y=abunds/1.197))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/1.197, ymax=(abunds+SE)/1.197))+
  scale_x_continuous("Canopy cover (%)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=14),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=14))




# Plot cat abundance by human population density (camera)

plotDensCam <- ggplot(data=cAbund, mapping=aes(x=hDensity, y=abunds/1.197))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/1.197, ymax=(abunds+SE)/1.197))+
  scale_x_continuous("Human population density (UNITS?)", 
                     breaks=c(0,10,20,30,40,50,60,70,80,90,100),
                     labels=c(0,10,20,30,40,50,60,70,80,90,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=14),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=14))



# Plot cat abundance by high school education (camera)

plotEduCam <- ggplot(data=cAbund, mapping=aes(x=eduHS*100,y=abunds/1.197))+
  geom_point(size=2)+
  geom_errorbar(mapping=aes(ymin=(abunds-SE)/1.197, ymax=(abunds+SE)/1.197))+
  scale_x_continuous("High school degree (%)", 
                     breaks=c(80,85,90,95,100),
                     labels=c(80,85,90,95,100))+
  scale_y_continuous("Scaled cat abundance",
                     breaks = c(0,.25,.5,.75,1.0),
                     labels=c(0,0.25,0.5,0.75,1.0))+
  coord_cartesian(ylim=c(0,1.6))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.title.y=element_text(margin=margin(0,20,0,0),size=10),
        axis.title.x=element_text(margin=margin(20,0,0,0),size=10))



# ==================================================================================*
# -------------------------------- Other plots -------------------------------------
# ==================================================================================*


# ----------------------------------------------------------------------------------*
# ---- Plots of urbanization by site ----
# ----------------------------------------------------------------------------------*


# Plot impervious surface % at study sites

plotStudyImp <- ggplot(data=covs)+
  geom_histogram(
    mapping=aes(x=imp, y=..count../sum(..count..)),
    binwidth=10,
    center=5
    )+
  scale_x_continuous(
    "Impervious surface (%)", 
    breaks=c(0,10,20,30,40,50,60,70,80,90,100),
    labels=c(0,10,20,30,40,50,60,70,80,90,100),
    limits=c(0,100)
    )+
  scale_y_continuous(
    "Proportion of sites",
    breaks=c(0,0.1,0.2,0.3,0.4,0.5),
    labels=c(0,0.1,0.2,0.3,0.4,0.5)
    )+
  coord_cartesian(ylim=c(0,0.5))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0),size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0),size=10))


# Plot impervious surface % at all Neighborhood Nestwatch sites

plotNNImp <- ggplot(data=NNSites)+
  geom_histogram(
    mapping=aes(x=imp, y=..count../sum(..count..)),
    binwidth=10,
    center=5
    )+
  scale_x_continuous(
    "Impervious surface (%)", 
    breaks=c(0,10,20,30,40,50,60,70,80,90,100),
    labels=c(0,10,20,30,40,50,60,70,80,90,100),
    limits=c(0,100)
    )+
  scale_y_continuous(
    "Proportion of sites",
    breaks=c(0,0.1,0.2,0.3,0.4,0.5),
    labels=c(0,0.1,0.2,0.3,0.4,0.5),
    limits=c(0,0.5)
    )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0),size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0),size=10))


# ----------------------------------------------------------------------------------*
# ---- Cat activity plots ----
# ----------------------------------------------------------------------------------*


# Plot histogram of times that cats were detected on cameras

plotCamDetection <- ggplot(data=camDet)+
  geom_histogram(
    mapping=aes(x=time/60, y=..count../sum(..count..)), bins=15
    )+
  scale_x_continuous(
    "Hour of day", 
    breaks=c(0,4,8,12,16,20,24),
    labels=c(0,4,8,12,16,20,24),
    limits=c(0,24)
    )+
  scale_y_continuous(
    "Proportion of observations",
    limits=c(0,.15)
    )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0),size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0),size=10))


# Plot histogram of times that cats were detected on transects

plotTransCats <- ggplot(data=transCat)+
  geom_histogram(
    mapping=aes(x=time/60, y=..count../sum(..count..)), bins=15
  )+
  scale_x_continuous(
    "Hour of day", 
    breaks=c(0,4,8,12,16,20,24),
    labels=c(0,4,8,12,16,20,24),
    limits=c(0,24)
  )+
  scale_y_continuous(
    "Proportion of observations",
    limits=c(0,0.6)
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0),size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0),size=10))


# Plot histogram of times that transects were conducted

plotTransAll <- ggplot(data=transDet)+
  geom_histogram(
    mapping=aes(x=time/60, y=..count../sum(..count..)), bins=15
    )+
  scale_x_continuous(
    "Hour of day", 
    breaks=c(0,4,8,12,16,20,24),
    labels=c(0,4,8,12,16,20,24),
    limits=c(0,24)
    )+
  scale_y_continuous(
    "Proportion of observations",
    limits=c(0,.35)
    )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_text(margin=margin(0,10,0,0),size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0),size=10))

