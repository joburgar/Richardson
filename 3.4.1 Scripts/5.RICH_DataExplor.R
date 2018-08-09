##############################################################
# 5.RICH_DataExplor.R
# Created by Joanna Burgar, 27-Apr-2018
# Exploring covariate data for Richardson GLMMs
#############################################################

###--- Load packages
library(reshape2)	# for formatting data frames
library(dplyr)		# for applying functions to subsets of data frames
library(ggplot2)	# for data visualization
library(tidyr)		# for data formatting functions
library(lattice)  # for graphing
#library(zoo)      # for date conversion
library(Hmisc)    # for running correlations

###--- Load data
setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

load("camop.RData") # load camera effort data
load("db3.RData") # load species detection data

glimpse(camop)
glimpse(db3)

DayLookup <- read.csv("DayLookup.csv", row.names=1)
glimpse(DayLookup)

DayLookup$Yr_Month <- as.factor(format(as.Date(DayLookup$Datep), "%Y-%m"))
Year <- DayLookup$Year
YrMonth <- DayLookup$YrMonth
YrWeek <- DayLookup$YrWeek
Day <- DayLookup$StudyDay


camloc <- read.csv("Richardson_CT_StationTable.csv")
glimpse(camloc)

colnames(camloc)[1] <- "Station"

db3$Line <- camloc$Line[match(db3$Station, camloc$Station)]
head(db3)

table(db3$Line, db3$Species)
sum(db3$count)

#############################################################
###--- aggregate species detection data by month and then week
# add year month and year week to db3
db3$Yr_Month <- as.factor(format(as.Date(db3$DateTimeOriginalp), "%Y-%m"))
db3$Yr_Week <- as.factor(format(as.Date(db3$DateTimeOriginalp), "%Y-%W"))
glimpse(db3)
#############################################################
###--- create monthly detection frame with Yr_Month the same as in db3
CA.month <- as.data.frame(CA.month)
CA.month$Location <- rownames(CA.month)
names(CA.month)

head(CA.month2)
CA.month2 <- gather(CA.month, "Yr_Month","Caribou",1:49) # create rows for each camera location by Yr_Month
glimpse(CA.month2)
CA.month2$Year <- substr(CA.month2$Yr_Month,2,5)
CA.month2$Month <- substr(CA.month2$Yr_Month,7,9)
CA.month2$YrMonth <- paste(CA.month2$Year,CA.month2$Month, sep="-")
CA.month2$Yr_Monthp <- as.yearmon(CA.month2$YrMonth, "%Y-%b")

CA.month2$Yr_Month <- as.factor(format(as.Date(CA.month2$Yr_Monthp), "%Y-%m"))
CA.month2$Year <- as.factor(CA.month2$Year)

levels(CA.month2$Yr_Month) # same levels as in db3
levels(db3$Yr_Month) # same levels as in CA.month2

CA.month2 <- CA.month2[order(CA.month2$Location),]
glimpse(CA.month2)
nrow(CA.month2) #3185

# create a data frame for GLMMs
CA.df <- CA.month2[c(1:5)] # keep only Location, Yr_Month, Year, Month and Caribou detection columns
CA.df <- CA.df %>%         # remove the survey start and end months (incomplete data)
  filter(Yr_Month!="2013-07") %>%
  filter(Yr_Month!="2017-07") 
nrow(CA.df)

###--- Create monthly detection dataframe for covariates
glimpse(db3)
db4 <- db3 %>%
  group_by(Location, Species, Yr_Month) %>% 
  summarise(sum(count, na.rm = TRUE))
colnames(db4) <- c("Location","Species","Yr_Month","Count")
table(db4$Yr_Month,db4$Species)
head(db4)

# create monthly counts for covariates for GLMMs
db4 <- as.data.frame(db4)
db4.month <- spread(db4, Species, Count)
db4.month[is.na(db4.month)] <- 0

# create combined covariates
db4.month$predator <- db4.month$`Grizzly Bear` + db4.month$Wolf # montly predator count column
db4.month$traffic <-  db4.month[,6] + db4.month[,7] + db4.month[,8] + db4.month[,9]  # monthly vehicle count column
head(db4.month)
summary(db4.month$traffic)

###--- Add covariates to caribou dataframe for GLMMs - for predators and traffic from detection data
CA.df <- left_join(CA.df, db4.month, by=c("Location","Yr_Month"))
glimpse(CA.df)
head(CA.df)
CA.df <- CA.df[,-6] # remove caribou from species detections (as CA caribou detections include NA for camera effort)
ncol(CA.df)
CA.df[,6:14][is.na(CA.df[,6:14])] <- 0 #change NA to 0 for all covariates (artefact of bringing over from detection data)

# add site level covariates 
CA.df <- left_join(CA.df, camdf[,c("Location","Park","HabitatType","Distance.to.HPAR","Elevation")], by="Location")
colnames(CA.df) <- c("Location","Yr_Month","Caribou","Year","Month","GB","MO","V_L","V_M","V_S","V_XL","WO","predator","traffic",
                     "park","habitat","distHPAR","elev")

CA.df$forest <- as.factor(ifelse(CA.df$habitat=="Open","Open","Forest"))
CA.df$Felev <- as.factor(ifelse(CA.df$elev<1200,1,
                                ifelse(CA.df$elev<1400,2,
                                       ifelse(CA.df$elev<1600,3,4))))
plot(CA.df$Felev,CA.df$elev)


CA.df$area <- as.factor(substr(CA.df$Location,1,6)) # linear area of grouped camera location distances (9 areas)
CA.df$site <- as.factor(substr(CA.df$Location,1,7)) # grid cell, camera may have moved locations during survey (61 sites)

# create covariate for west, east or on road
names(camdf)
road.utm <- camdf[camdf$Distance.to.HPAR<11,c(1,3)] # first create vector of just road camera easting coordinates
road.utm$area <- substr(road.utm$Location,1,6) # linear area of grouped camera location distances (9 areas)

CA.df$road.utmx <- road.utm$utm_x[match(CA.df$area, road.utm$area)]
CA.df$utmx <- camdf$utm_x[match(CA.df$Location, camdf$Location)]
CA.df$road <- as.factor(ifelse(CA.df$distHPAR>10 & CA.df$utmx<=CA.df$road.utmx, "East",
                               ifelse(CA.df$distHPAR>10 & CA.df$utmx>=CA.df$road.utmx, "West","Road")))
CA.df$road <- relevel(CA.df$road,"Road") # relevel factor with "Road" as reference level

plot(CA.df$road,CA.df$traffic)

# create traffic volume to be used with distance to HPAR - use traffic volume at road camera for entire linear area
head(CA.df)
traffic.vol <- CA.df[CA.df$road=="Road",c("Yr_Month","traffic","area","road")]
glimpse(traffic.vol)
head(traffic.vol)
CA.df$trf.vol <- traffic.vol$traffic[match(CA.df$area, traffic.vol$area)]


# add temporal covariates
CA.df$Season <- as.factor(ifelse(CA.df$Month=="Mar"|CA.df$Month=="Apr"|CA.df$Month=="May","spring",
                                 ifelse(CA.df$Month=="Jun"|CA.df$Month=="Jul"|CA.df$Month=="Aug","summer",
                                        ifelse(CA.df$Month=="Sep"|CA.df$Month=="Oct"|CA.df$Month=="Nov","fall","winter"))))
CA.df$Season <- relevel(CA.df$Season,"summer") # relevel factor with "summer" as reference level

head(CA.df)
summary(CA.df)

# change columns with class(double) to class(integer)
cols_to_change = c(6:14)
for(i in cols_to_change){
  class(CA.df[, i]) = "integer"
}

glimpse(CA.df)

###############################################################################################

# find correlated variables
names(camdf)
camdf.cor <- as.matrix(camdf[c(31:33)])
rcorr(camdf.cor, type="pearson")
# distance to HPAR correlated with distance to LNR (0.86) and elevation (0.77)

names(CA.df)
CA.df.cor <- as.matrix(CA.df[c(7,29,31)])
rcorr(CA.df.cor, type="pearson")

#############################################################
### Basic data exploration 

names(CA.df)

#Outliers
par(mfrow = c(1, 2))
boxplot(CA.df$Caribou, 
        main = "Monthly Caribou Detections",
        horizontal=TRUE)

boxplot(CA.df$predator, 
        main = "Monthly Pedator Detections",
        horizontal=TRUE)

par(mfrow = c(1, 5), mar = c(4, 3, 3, 2))
dotchart(CA.df$Caribou, main = "Monthly Caribou Detections")
dotchart(CA.df$traffic, main = "Monthly Traffic Volume")
dotchart(CA.df$predator, main = "Monthly Predator Detections")
dotchart(CA.df$distHPAR, main = "Distance to HPAR")
dotchart(CA.df$elev, main = "Elevation")

# frequency plots to see how caribou detections vary by site/area
dotplot(CA.df$site,db4.month$Caribou)
dotplot(CA.df$area,db4.month$Caribou)
par(mfrow = c(1, 2))
plot(as.factor(CA.df$road), CA.df$Caribou)
plot(as.factor(CA.df$Felev), CA.df$Caribou) # where 1=1000-1200 and each level is 200 m increment


###--- basic relationship plots
plot(CA.df$elev,CA.df$distHPAR)


par(mfrow = c(2,3), mar = c(4, 3, 3, 2))
# landscape covariates
plot(CA.df$distHPAR,CA.df$Caribou)
plot(CA.df$elev,CA.df$Caribou)
plot(CA.df$habitat,CA.df$Caribou)
plot(as.factor(CA.df$forest),CA.df$Caribou)
plot(CA.df$park,CA.df$Caribou)
plot(as.factor(CA.df$road), CA.df$Caribou)

# temporal covariates
par(mfrow = c(1,2), mar = c(4, 3, 3, 2))
plot(as.factor(CA.df$Year),CA.df$Caribou)
plot(as.factor(CA.df$Season),CA.df$Caribou)

# behavioural covarariates
par(mfrow = c(1,2), mar = c(4, 3, 3, 2))
plot(CA.df$predator,CA.df$Caribou)
plot(CA.df$traffic,CA.df$Caribou)

##############################################
###--- Plotting possible interactions
##-- Caribou ~ Distance to HPAR
xyplot(Caribou ~  distHPAR| factor(park),
       data = CA.df, 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         panel.loess(x,y,col=1,lwd=2) #Add smoother
         #panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  distHPAR| factor(Season),
       data = CA.df, 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         panel.loess(x,y,col=1,lwd=2) #Add smoother
        #panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  distHPAR| factor(Season),
       data = CA.df[CA.df$park=="NaNPR",], 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  distHPAR| factor(Season),
       data = CA.df[CA.df$park=="NNPR",], 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  distHPAR| factor(forest),
       data = CA.df, 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         panel.loess(x,y,col=1,lwd=2) #Add smoother
         #panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  distHPAR| factor(road),
       data = CA.df, 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })

xyplot(Caribou ~  distHPAR| factor(Felev),
       data = CA.df, 
       xlab = "Distance to HPAR (m)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })

##-- Caribou ~ Predator Detections
xyplot(Caribou ~  predator| factor(park),
       data = CA.df, 
       xlab = "Monthly Predator Detections",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  predator| factor(Season),
       data = CA.df, 
       xlab = "Monthly Predator Detections",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  predator| factor(forest),
       data = CA.df, 
       xlab = "Monthly Predator Detections",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  predator| factor(road),
       data = CA.df, 
       xlab = "Monthly Predator Detections",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  predator| factor(Felev),
       data = CA.df, 
       xlab = "Monthly Predator Detections",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


##############################################
##-- Caribou ~ Traffic Volume (using only road cameras)
CA.road <- droplevels(subset(CA.df, road == "Road"))

xyplot(Caribou ~  traffic| factor(park),
       data = CA.road, 
       xlab = "Traffic Volume (Monthly Vehicle Passes)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  traffic| factor(Season),
       data = CA.road, 
       xlab = "Traffic Volume (Monthly Vehicle Passes)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  traffic| factor(forest),
       data = CA.road, 
       xlab = "Traffic Volume (Monthly Vehicle Passes)",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


###--- Exploring traffic volume in each park by vehicle type
xyplot(Caribou ~  V_XL | factor(park),
       data = CA.road,
       xlab = "XL Vehicle Volume",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })

xyplot(Caribou ~  V_L | factor(park),
       data = CA.road,
       xlab = "Large Vehicle Volume",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })


xyplot(Caribou ~  V_M | factor(park),
       data = CA.road,
       xlab = "Medium Vehicle Volume",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })

xyplot(Caribou ~  V_S | factor(park),
       data = CA.road,
       xlab = "Small Vehicle Volume",
       ylab = "Monthly Caribou Detections",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       scales = list(alternating = TRUE, 
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel=function(x,y){
         panel.grid(h=-1, v= 2)
         panel.points(x, y, col = 1)
         #panel.loess(x,y,col=1,lwd=2) #Add smoother
         panel.abline(lm(y~x))        #Add regression line
       })