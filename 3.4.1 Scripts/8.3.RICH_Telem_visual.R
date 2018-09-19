#############################################################
# 8.3.RICH_Telem_visual.R
# created by Joanna Burgar, 19-Sept-2018
# script for visualizing telemetry of collared caribou
#############################################################

###--- Code to take telemetry data from csv text file (i.e., downloaded from ArcGIS and animate)

###--- Load packages
library(plyr) #  might not be necessary
library(dplyr)    # for viewing and manipulating data
#library(tidyr)    # for manipulating data  - might not be necessary
#library(reshape2) # for manipulating data  - might not be necessary
library(anipaths) # for making animated graphics
library(ggmap)    # for background of maps

# Load data
setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

ca.df <- read.table("RICH_caribouGPS_SAoverlap13Nov-16Jan.txt", header=TRUE, sep=",") 
head(ca.df)
glimpse(ca.df)

###--- Add columns for the Date/Time in POSIX format:
ca.df$Datep <- as.POSIXct(strptime(ca.df$Date, format = "%d/%m/%Y"))
ca.df$DateTime <- paste(ca.df$Datep, ca.df$Time)
ca.df$DateTimep <- as.POSIXct(strptime(ca.df$DateTime, format = "%Y-%m-%d %H:%M:%S"))

ca.df$Indiv <- as.factor(ca.df$Individual)
summary(ca.df)

ca.df1 <- ca.df %>% 
  filter(Individual="2228")

min(ca.df$Xcoord)
max(ca.df$Ycoord)

min(ca.df$LongWGS84)
max(ca.df$LatWGS84)

summary(ca.df2)
background.latlong <- list(location = c(min(ca.df$LongWGS84), min(ca.df$LatWGS84)), 
                       zoom = "auto", 
                       maptype = "satellite")


animate_paths(paths = ca.df, 
              delta.t = "day",
              coord = c("LongWGS84","LatWGS84"), # set up as long and then lat
              Time.name = "DateTimep",
              covariate.colors = RColorBrewer::brewer.pal(n = 2, "RdYlGn"),
              ID.name = "Indiv",
              background = background.latlong)
             
#file:///C:/Users/JBurgar/Google%20Drive/Richardson%20Wildfire%20Project/3.%20Data/3.4%20Data%20Analysis/3.4.2%20Input/index.html
