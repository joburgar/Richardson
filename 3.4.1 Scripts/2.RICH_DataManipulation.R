##############################################################
# 2.RICH_DataManipulation.R
# Created by Joanna Burgar, 26-Apr-2018
# Formatting Richardson species detection data to use in analyses
#############################################################

###--- Load packages
library(plyr)
library(dplyr)		# for applying functions to subsets of data frames
library(camtrapR) # for creating species matrices

###--- Import csv file - output record table from Camelot, spreadsheet of camera operability
setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory
load("camop.RData") # load camdf data

# Species detection data - load record table from Camelot
# exported 3 record tables from Camelot - 10, 30 and 60 min thresholds - identical for independent detection events

db1 <- read.csv("record-table_2018-04-26_Richardson_30min.csv") 
glimpse(db1) # check for correct loading of data - to see data classes and initial values
summary(db1) # overview of data and check for NAs


#############################################################
###--- Format date-time for R

###--- Add columns for the Date/Time in POSIX format:
db1$DateTimeOriginalp <- as.POSIXct(strptime(db1$DateTimeOriginal, format = "%Y-%m-%d %H:%M:%S"))
db1$Datep <- as.POSIXct(strptime(db1$Date, format = "%Y-%m-%d"))
db1$Timep <- as.POSIXct(strptime(db1$Time, format = "%H:%M:%S"))
db1$Timep <- strftime(db1$Timep, format="%H:%M:%S")

###--- Create Year, Month, Week columns (factors)
db1$Year <- as.factor(format(as.Date(db1$DateTimeOriginalp),"%Y"))
db1$Month <- as.factor(months(db1$Datep, abbreviate=TRUE))
db1$Week <- as.factor(format(as.Date(db1$Datep),"%V"))

glimpse(db1)
summary(db1)
head(db1)
names(db1)

###--- Subset data to columns used in analyses
db2 <- db1[c(1,4,15:20)]
head(db2)
summary(db2)
names(db2)

db2 <- db2[order(db2$DateTimeOriginalp),]
tail(db2) #for some reason one entry for Lepus americanus is coming up with NA for DateTimep and Year - manually fixing for the moment

db2[751,3] <- "2018-03-11 02:31:08"
db2[751,6] <- 2018

#############################################################
###--- Create dataframe for only species of interest
# check for errors in "Species" variable

levels(db2$Species)

# subset data for known mammal species (other than humans)
db3 <- db2 %>%
  filter(Species!="Bird spp."& 
           Species!="Unknown species"&
           Species !="No Animal"&
           Species != "Homo sapiens")

summary(db2$Species)
summary(db3$Species)

glimpse(db3)


# check total number of records per species and for each year
db3$count <- 1 #create a count of 1 for each record

db3 %>% group_by(Species) %>% summarise(sum(count)) # table with # records per "Species"
#Species                 `sum(count)`
#1 Alces alces                       8.
#2 Canis lupus                       5.
#3 Lepus americanus                153.
#4 Lontra canadensis                 7.
#5 Lynx canadensis                  25.
#6 Martes americana                 14.
#7 Martes pennanti                   3.
#8 Odocoileus virginianus            1.
#9 Rangifer tarandus                29.
#10 Tamiasciurus hudsonicus           9.
#11 Vulpes vulpes                     3

# 11 species detected, range from 1 detection for O. virginianus to 153 detections of L. americanus

db3 %>% group_by(Year, Species) %>% summarise(sum(count)) # table with # records per "Species" and Year

#############################################################
###--- Housekeeping - cleaning up files
###--- Saving only db3 dataframe
rm(camdf, camop, db1, db2)
save.image(file="db3.RData") 
