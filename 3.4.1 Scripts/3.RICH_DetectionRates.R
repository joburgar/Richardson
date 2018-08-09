##############################################################
# 3. RICH_DetectionRates.R
# Created by Joanna Burgar, 26-Apr-2018
# Creating Richardson species rates & species matrices
#############################################################

###--- Load packages
library(reshape2)
library(dplyr)

setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

load("camop.RData") # load camera effort data
names(camdf)

load("db3.RData") # load species detection data
names(db3)
View(db3)

#############################################################
###--- not all cameras were operational for entire study period
###--- range of detection timing depend on location

# convert DateStart (but note that this doesn't have time set, so will treat as midnight)
db3$DateStart <- camdf$setup_date[match(db3$Station, camdf$SiteID)] 
db3$DateStart <- (strptime(db3$DateStart, "%Y-%m-%d", tz="MST"))
db3$DateStart <- as.POSIXct(db3$DateStart)
DateStart <- min(db3$DateStart) #  "2017-11-11 MST"
summary(db3)

# calculate a unique day for each day of study
# taking straight difference will include partial days, so is really 24-hour periods from time first camera set
# using "floor" so it doesn't round up to next day
db3$StudyDay <- floor(as.numeric(difftime(db3$Datep,min(db3$DateStart),units="days")))
db3$StudyDay <- db3$StudyDay+1
max(db3$StudyDay) # 151 study days

DayLookup <- array(0,dim=c(max(db3$StudyDay),2))
colnames(DayLookup) <- c("StudyDay", "Datep")
dim(DayLookup)
max(db3$DateTimeOriginalp) - DateStart
  
DayLookup <- as.data.frame(DayLookup)
DayLookup$StudyDay <- seq.int(max(db3$StudyDay))
DayLookup$Datep <- seq(DateStart, max(db3$DateTimeOriginalp), by="days")

DayLookup$Datep <- as.Date(DayLookup$Datep, format="%Y-%m-%d")

DayLookup$Year <- as.factor(format(as.Date(DayLookup$Datep),"%Y"))
DayLookup$Month <-  as.factor(months(DayLookup$Datep, abbreviate=TRUE))
DayLookup$Week <- as.factor(format(as.Date(DayLookup$Datep),"%V"))

summary(DayLookup)
glimpse(DayLookup)

write.csv(DayLookup,"DayLookup.csv") # to be used as input data in other analyses

unique(db3$Station) # 29 sites with mammal data
summary(camdf$SiteID) # 30 CT sites
dim(camop) # 30 by 152

#############################################################
###--- function to create species detection matrix

fn.sp.dtnmtrx <- function(D=D, Species=Species, effort=effort){
 
  D.ord <- D[order(D$Station),]
  Nsp <- D[which(D$Species==Species),]
  with(Nsp, table(StudyDay, count))
  with(Nsp, table(Station, count))
  
  length(unique(Nsp$StudyDay)) 
  sum(Nsp$count) 
  
  Sp.tmp <- dcast(Nsp,Station~StudyDay,fun.aggregate = sum, value.var="count",fill=0) 
  dim(Sp.tmp) 
  
  row.names(Sp.tmp) <- Sp.tmp$Station
  Sp.tmp <- Sp.tmp[,-1]
  
  study.days <- 1:max(D$StudyDay)
  Sp.missing <- study.days[-which(study.days %in% names(Sp.tmp))]
  Sp.zeros <- as.data.frame.matrix(matrix(0,nrow=nrow(camop),ncol=length(Sp.missing)))
  names(Sp.zeros)<- as.character(Sp.missing) 
  row.names(Sp.zeros) <- row.names(camop)
  
  Sp.tmp0 <- as.data.frame(matrix(0,ncol = length(Sp.tmp), nrow = nrow(Sp.zeros) - nrow(Sp.tmp)))
  names(Sp.tmp0) <- names(Sp.tmp)
  row.names(Sp.tmp0) <- row.names(Sp.zeros)[-which(row.names(Sp.zeros) %in% row.names(Sp.tmp))]
  Sp.tmp2 <- rbind(Sp.tmp,Sp.tmp0)
  Sp.tmp2 <- Sp.tmp2[order(as.numeric(row.names(Sp.tmp2))),]

  Sptmp <- data.frame(Sp.tmp2,Sp.zeros,check.names=FALSE)
  Sp.mat <- Sptmp[,order(as.numeric(names(Sptmp)))]
  Sp.mat2 <- Sp.mat[order(rownames(Sp.mat)),]
  Sp.mat2 <- as.data.frame(Sp.mat2)
  SpMatrix <- Sp.mat2
  
  effort <- effort[order(rownames(effort)), ] 
  
  for (SITE in 1:nrow(SpMatrix)) {
    for (DAY in 1:ncol(SpMatrix)) {
      SpMatrix[SITE,DAY] <- ifelse(effort[SITE,DAY] == 0, NA, SpMatrix[SITE,DAY])
    } }
  
  return(SpMatrix=SpMatrix)
  }

#############################################################

levels(db3$Species)
summary(db3$Species)

###--- species detection matrices
# disregard warning messages re:NAs introduced by coercion

# Caribou
CA.matrix <- fn.sp.dtnmtrx(D=db3, Species="Rangifer tarandus", effort=camop)
dim(CA.matrix) # 30 by 151 - check 
sum(CA.matrix, na.rm = TRUE) #29 - same as in db3
summary(db3$Species[db3$Species=="Rangifer tarandus"]) #29 in db3

# Lynx
LY.matrix <- fn.sp.dtnmtrx(D=db3, Species="Lynx canadensis", effort=camop)
dim(LY.matrix) # 30 by 151 - check 
sum(LY.matrix, na.rm = TRUE) #25 - same as in db3
summary(db3$Species[db3$Species=="Lynx canadensis"]) #25 in db3

# Marten
MA.matrix <- fn.sp.dtnmtrx(D=db3, Species="Martes americana", effort=camop)
dim(MA.matrix) # 30 by 151 - check 
sum(MA.matrix, na.rm = TRUE) #14 - same as in db3
summary(db3$Species[db3$Species=="Martes americana"]) #14 in db3

# Moose
MO.matrix <- fn.sp.dtnmtrx(D=db3, Species="Alces alces", effort=camop)
dim(MO.matrix) # 30 by 151 - check 
sum(MO.matrix, na.rm = TRUE) #8 - same as in db3
summary(db3$Species[db3$Species=="Alces alces"]) #8 in db3

# Snowshoe Hare
SH.matrix <- fn.sp.dtnmtrx(D=db3, Species="Lepus americanus", effort=camop)
dim(SH.matrix) # 30 by 151 - check 
sum(SH.matrix, na.rm = TRUE) #153 - same as in db3
summary(db3$Species[db3$Species=="Lepus americanus"]) #153 in db3

# Wolf
WO.matrix <- fn.sp.dtnmtrx(D=db3, Species="Canis lupus", effort=camop)
dim(WO.matrix) # 30 by 151 - check 
sum(WO.matrix, na.rm = TRUE) #5 - same as in db3
summary(db3$Species[db3$Species=="Canis lupus"]) #5 in db3


#############################################################
#############################################################
###--- Housekeeping - cleaning up files
###--- Saving only daylookup and species matrices - caribou, lynx, marten, moose, snowshoe hare and wolf
ls()

rm(list= ls()[!(ls() %in% c('CA.matrix','LY.matrix','MA.matrix','MO.matrix','SH.matrix','WO.matrix',
                            'DayLookup'))])
save.image(file="SpMatrix.RData") 
#load("SpMatrix.RData")



#############################################################
### detector operation schedule in the form of an image plot
library(scrbook)
oper <- camop
oper[is.na(oper)] <- 0
# dark blue is not operational and light blue is operational
image(1:ncol(oper),1:nrow(oper),rot(oper), xlab="Study Day", ylab="Number of Detectors",col=topo.colors(2))

###--- Caribou detections
dim(CA.matrix)

# some summaries on caribou daily detections
rowSums(CA.matrix,na.rm=T)
sum(CA.matrix,na.rm=T) #29
dotchart(rowSums(CA.matrix,na.rm=T),pch=16,ylab="Camera Station",xlab="Total Caribou Count")
barplot(sort(rowSums(CA.matrix,na.rm=T)))
plot(colSums(CA.matrix[,1:151],na.rm=T), xlab="Study Day", ylab="Caribou Detections")

DayLookup[which(DayLookup$StudyDay==(1)),] # 2017-11-11 - 1st study day
DayLookup[which(DayLookup$StudyDay==(151)),] # 2018-04-10 - last study day

# Cameras put out 11-13 Nov so start "survey period" 14 Nov
# Blue lines = Nov 14 - Jan 17 for ~2 months (64 days) and maximise detections

# Red lines = post-calving (1 June - 31 Aug)
# Green lines = rut and immediately post rut (15 Sept - 15 Dec)

DayLookup[which(DayLookup$Datep=="2017-11-14"),] # 4
DayLookup[which(DayLookup$StudyDay==(68)),] # 13 Jan
abline(v=c(4,68),col="blue")
