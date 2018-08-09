##############################################################
# 4. RICH_Spatial-Temporal.R
# Created by Joanna Burgar, 26-Apr-2018
# Spatial and temporal patterns for species use
#############################################################

###--- Load packages
#library(reshape2)	# for formatting data frames
library(plyr)
library(dplyr)		# for data manipulation
library(camtrapR) # for visualizing species detection data

###--- Load data
setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

load("camop.RData") # load camera effort data
load("db3.RData") # load species detection data

setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.3 Output") # set working directory

#############################################################
###--- Visualization of species detection data
# uses camtrapR and dependent packages
names(camdf)
colnames(camdf)[1] <- "Station"
names(db3)

###--- Spatial visualization
sp.det.map <- detectionMaps(CTtable = camdf,
                            recordTable = db3,
                            Xcol = "utm_x",
                            Ycol = "utm_y",
                            stationCol = "Station",
                            speciesCol = "Species",
                            writePNG = TRUE,
                            plotDirectory = getwd(),
                            plotR = FALSE,
                            printLabels = FALSE,
                            richnessPlot = FALSE,
                            addLegend = TRUE)


###--- Temporal visualization

# Activity Histograms for each focal species
act.hist <- activityHistogram(recordTable = db3,
                                   allSpecies = TRUE,
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTimeOriginalp",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   plotR = TRUE,
                                   writePNG = FALSE) 


# Activity density graphs for each focal species
act.dens <- activityDensity(recordTable = db3,
                                 allSpecies = TRUE,
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTimeOriginalp",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 plotR = TRUE,
                                 writePNG = FALSE) 


# Activity raidal graphs or each focal species
act.rad <- activityRadial(recordTable = db3,
                               allSpecies = TRUE,
                               speciesCol = "Species",
                               recordDateTimeCol = "DateTimeOriginalp",
                               recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                               byNumber = TRUE,
                               plotR = TRUE,
                               writePNG = FALSE) 



#############################################################
###--- Visualization of overlapping species detection data
summary(db3$Species)

# Caribou and Wolf
overlap.CA.WO <- activityOverlap(recordTable = db3,
                                 speciesA = "Rangifer tarandus",
                                 speciesB = "Canis lupus",
                                 allSpecies = TRUE,
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTimeOriginalp",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 plotR = TRUE,
                                 writePNG = TRUE,
                                 plotDirectory = getwd(),
                                 addLegend = TRUE,
                                 legendPosition = "topleft") 
# overlap Dhat=0.56

overlap.LA.LY <- activityOverlap(recordTable = db3,
                                 speciesA = "Lepus americanus",
                                 speciesB = "Lynx canadensis",
                                 allSpecies = TRUE,
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTimeOriginalp",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 plotR = TRUE,
                                 writePNG = TRUE,
                                 plotDirectory = getwd(),
                                 addLegend = TRUE,
                                 legendPosition = "topleft") 
# overlap Dhat=0.48


###############################################################################
# Caribou and Moose
overlap.CA.MO <- activityOverlap(recordTable = db3,
                                 speciesA = "Rangifer tarandus",
                                 speciesB = "Alces alces",
                                 allSpecies = TRUE,
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTimeOriginalp",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 plotR = TRUE,
                                 writePNG = FALSE,
                                 #plotDirectory = getwd(),
                                 addLegend = TRUE,
                                 legendPosition = "topleft") 
# overlap Dhat=0.42


# Moose and Wolf
overlap.MO.WO <- activityOverlap(recordTable = db3,
                                 speciesA = "Alces alces",
                                 speciesB = "Canis lupus",
                                 allSpecies = TRUE,
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTimeOriginalp",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 plotR = TRUE,
                                 writePNG = TRUE,
                                 plotDirectory = getwd(),
                                 addLegend = TRUE,
                                 legendPosition = "topleft") 
# overlap Dhat=0.41
