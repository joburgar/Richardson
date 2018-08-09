##############################################################
# 1. RICH_CameraData.R
# Created by Joanna Burgar, 26-Apr-2018
# Formatting Richardson camera data to use in analyses
#############################################################

###--- Load packages
library(plyr)
library(dplyr)    # for viewing and manipulating data
library(camtrapR) # for creating camera location operability matrix

###--- Import csv file - output record table from Camelot, spreadsheet of camera operability
setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

camloc <- read.csv("Richardson_CT_StationTable.csv")

glimpse(camloc) # check for correct loading of data - to see data classes and initial values
summary(camloc) # overview of data and check for NAs

levels(camloc$SiteID) # 30 sites
nrow(camloc)
summary(camloc) # 18 sites on game trails, 12 on seismic lines - all outside of burns

###--- Add columns for the Dates in POSIX format:
names(camloc)
camloc$setup_date <- as.POSIXct(strptime(camloc$setup_date, format = "%d-%b-%y"))
camloc$retrieval_date <- as.POSIXct(strptime(camloc$retrieval_date, format = "%d-%b-%y"))

camloc$Problem1_from <- as.POSIXct(strptime(camloc$Problem1_from, format = "%d-%b-%y"))
camloc$Problem1_to <- as.POSIXct(strptime(camloc$Problem1_to, format = "%d-%b-%y"))
camloc$Problem2_from <- as.POSIXct(strptime(camloc$Problem2_from, format = "%d-%b-%y"))
camloc$Problem2_to <- as.POSIXct(strptime(camloc$Problem2_to, format = "%d-%b-%y"))
camloc$Problem3_from <- as.POSIXct(strptime(camloc$Problem3_from, format = "%d-%b-%y"))
camloc$Problem3_to <- as.POSIXct(strptime(camloc$Problem3_to, format = "%d-%b-%y"))

glimpse(camloc)

camdf <-camloc[,c(1:5,9:14)] # remove covariates and extra columns
nrow(camdf) # 30 cameras
head(camdf)
View(camdf)

mindate <- aggregate(setup_date~SiteID, camdf, function(x) min(x))
maxdate <- aggregate(retrieval_date~SiteID, camdf, function(x) max(x))

mindate <- mindate[order(mindate$SiteID),] # order by Location
mindate
maxdate <- maxdate[order(maxdate$SiteID),] # order by Location
maxdate

min(maxdate$retrieval_date - mindate$setup_date) # 148 days
max(maxdate$retrieval_date - mindate$setup_date) # 150 days
mean(maxdate$retrieval_date - mindate$setup_date) # 149 days
# range of 148-150 days for cameras between setup and retrieval

camdf[order(camdf$Problem1_from),]
summary(camdf)

###--- Camera trap location operability matrix
camop <- cameraOperation(camdf, 
                stationCol = "SiteID", 
                setupCol = "setup_date", 
                retrievalCol = "retrieval_date", 
                hasProblems = TRUE,
                dateFormat = "%Y-%m-%d", 
                writecsv = TRUE,
                outDir = getwd())


View(camop)
#Legend: NA: camera(s) not set up, 0: camera(s) not operational, 1 (or higher): number of operational camera(s)
dim(camop) # 30 locations x 152 days

camop[camop == 0] <- NA # changes all 0 to NA, as 0 is not operational

camop <- camop[order(rownames(camop)), ]  # order alphabetically by Location (same as camdf)
camdf <- camdf[order(camdf$SiteID), ]   # order alphabetically by Location (same as camop)

sum(camop, na.rm=TRUE) #4231 camera trap days - total effort for survey to date

#############################################################
###--- Housekeeping - cleaning up files
###--- Saving only camop matrix and camdf covariate dataframe

rm(list= ls()[!(ls() %in% c('camop','camdf'))])
save.image(file="camop.RData") 
