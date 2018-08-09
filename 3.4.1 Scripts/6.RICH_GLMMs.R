##############################################################
# HPAR_GLMMs.R
# Created by Joanna Burgar, 04-Mar-2018
# Creating HPAR GLMMs
#############################################################

###--- Load packages
library(glmmTMB) # for running GLMMs
citation("bbmle")
library(bbmle) #AICtab function

###--- Can skip lines of code (21-155) if already run HPAR_DataExplor.R ---####

library(reshape2)	# for formatting data frames
library(plyr)     # for manipulation data
library(dplyr)		# for applying functions to subsets of data frames
library(ggplot2)	# for data visualization
library(tidyr)		# for data formatting functions
library(lattice)  # for graphing
library(zoo)      # for date conversion
library(Cairo)    # for creating publication ready graphics

###--- Load data - 
setwd("C:/Users/JBurgar/Documents/R/Analysis/HPAR/") # set working directory

CA.month <- read.csv("CA.month.GLMM.csv", row.names=1) # load in monthly caribou detection matrix (with NA for month camera inactive)
#CA.week <- read.csv("CA.week.GLMM.csv", row.names=1) # load in weekly caribou detection matrix (with NA for week camera inactive)

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

#############################################################
###--- aggregate species detection data by month and then week
# add year month and year week to db3
db3$Yr_Month <- as.factor(format(as.Date(db3$DateTimep), "%Y-%m"))
db3$Yr_Week <- as.factor(format(as.Date(db3$DateTimep), "%Y-%W"))
glimpse(db3)
#############################################################
###--- create monthly detection frame with Yr_Month the same as in db3
CA.month <- as.data.frame(CA.month)
CA.month$Location <- rownames(CA.month)
names(CA.month)
ncol(CA.month)

CA.month2 <- gather(CA.month, "Yr_Month","Caribou",1:49) # create rows for each camera location by Yr_Month
CA.month2$Year <- substr(CA.month2$Yr_Month,2,5)
CA.month2$Month <- substr(CA.month2$Yr_Month,7,9)
CA.month2$YrMonth <- paste(CA.month2$Year,CA.month2$Month, sep="-")
CA.month2$Yr_Monthp <- as.yearmon(CA.month2$YrMonth, "%Y-%b")

CA.month2$Yr_Month <- as.factor(format(as.Date(CA.month2$Yr_Monthp), "%Y-%m"))

CA.month2$Year <- as.factor(strptime(CA.month2$Year, format = "%Y"))
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

summary(CA.df)

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

#############################################################
###--- standardize the integer data for regression
CA.df$distHPAR_scaled <- scale(CA.df$distHPAR, center=TRUE)
CA.df$distHPAR_km <- CA.df$distHPAR/1000 
summary(CA.df$distHPAR_km) # similar value range to Caribou, use km rather than scaled data

CA.df$predator_scaled <- scale(CA.df$predator, center=TRUE)
summary(CA.df$predator_scaled) # similar to predator value range, prefer use of unscaled
summary(CA.df$predator) # similar value range to Caribou, can use count data rather than scaled data
summary(CA.df$Caribou)

summary(CA.df)

CA.df$traffic_log <- log(CA.df$traffic+1)
summary(CA.df$traffic_log) # similar value range to Caribou

###--- GLMM poisson regressions
###--- with site and month as random factors
### Modelling with glmmTMB


hist(CA.df$Caribou, breaks=20)
m0 <- glm(Caribou ~ distHPAR_km, family=poisson,data=CA.df)

m1.1 <- glmmTMB(Caribou ~ distHPAR_km + (1|site), family=poisson, data=CA.df)
m1.2 <- glmmTMB(Caribou ~ distHPAR_km + (1|area), family=poisson, data=CA.df)
m1.3 <- glmmTMB(Caribou ~ distHPAR_km + (1|Month), family=poisson, data=CA.df)

m1.4 <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), family=poisson, data=CA.df)
m1.5 <- glmmTMB(Caribou ~ distHPAR_km + (1|area) + (1|Month), family=poisson, data=CA.df)

# use AIC to determine what random effects to include

AICtab(m1.5,m1.4,m1.3,m1.2,m1.1,m0, weights=TRUE) # model with random effects for both site and month is best


# to determine if poisson or negative binomial fits the data better
m1.pos <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), data=CA.df, family = poisson)
m1.posz <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), data=CA.df, zi=~1, family = poisson)

m1.nb <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), data=CA.df, family = nbinom2)
m1.nbz <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), data=CA.df, zi=~1, family = nbinom2)

anova(m1.nbz,m1.nb,m1.posz,m1.pos) # best model structure is negative binomial without zero inflation
AICtab(m1.nbz,m1.nb,m1.posz,m1.pos, weights=TRUE) 

###--- create regression models for information-theoretic approach
# null model
m0 <- glmmTMB(Caribou ~ 1 + (1|site) + (1|Month), data=CA.df, family = nbinom2)

###--- behavioural factors
mb.1 <- glmmTMB(Caribou ~ distHPAR_km + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.2 <- glmmTMB(Caribou ~ traffic_log + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.3 <- glmmTMB(Caribou ~ predator +  (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.4 <- glmmTMB(Caribou ~ distHPAR_km + predator + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.5 <- glmmTMB(Caribou ~ distHPAR_km + traffic_log + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.5b <- glmmTMB(Caribou ~ distHPAR_km*traffic_log + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.6 <- glmmTMB(Caribou ~ distHPAR_km + traffic_log + predator + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.7 <- glmmTMB(Caribou ~ MO +  (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.8 <- glmmTMB(Caribou ~ distHPAR_km + MO +  (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.9 <- glmmTMB(Caribou ~ distHPAR_km + MO + predator+ (1|site) + (1|Month), data=CA.df,family = nbinom2)
mb.10 <- glmmTMB(Caribou ~ distHPAR_km + MO + traffic_log +  (1|site) + (1|Month), data=CA.df,family = nbinom2)
anova(mb.3,mb.2,mb.1,m0) 

AICtab(mb.9,mb.6, mb.10,mb.8, mb.7,mb.5,mb.5b,mb.4,mb.3,mb.2,mb.1,m0, weights=TRUE) # models with interaction distance to road and traffic volume top model
summary(mb.5b)


###--- habitat/landscape factors
mh.1 <- glmmTMB(Caribou ~ park + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mh.2 <- glmmTMB(Caribou ~ forest + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mh.3 <- glmmTMB(Caribou ~ road + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mh.4 <- glmmTMB(Caribou ~ Felev + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mh.5 <- glmmTMB(Caribou ~ park + forest + road + Felev + (1|site) + (1|Month), data=CA.df,family = nbinom2)

anova(mh.5,mh.1,m0) # park no better than null
anova(mh.5,mh.2,m0) # forest no better than null
anova(mh.5,mh.3,m0) # road no better than null
anova(mh.5,mh.4,m0) # Felev no better than null

AICtab(mh.5,mh.4,mh.3,mh.2,mh.1,m0, weights=TRUE) # no model better than null 

###--- temporal factors (remove month random effect)
mt.1 <- glmmTMB(Caribou ~ Season + (1|site) + (1|Month),data=CA.df,family = nbinom2)
mt.2 <- glmmTMB(Caribou ~ as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mt.3 <- glmmTMB(Caribou ~ Season + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)

anova(mt.3,mt.2,mt.1) 
AICtab(mt.3,mt.2,mt.1, weights=TRUE) # Year (mt.3) is the best option (season and year)

# model with best of all three components
mfull.0 <- glmmTMB(Caribou ~ distHPAR_km*traffic_log +  Season + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mfull.1 <- glmmTMB(Caribou ~ distHPAR_km*Season +traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mfull.2 <- glmmTMB(Caribou ~ distHPAR_km*park + traffic_log + Season + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mfull.3 <- glmmTMB(Caribou ~ distHPAR_km*forest + traffic_log + Season + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mfull.4 <- glmmTMB(Caribou ~ distHPAR_km*road + traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
mfull.5 <- glmmTMB(Caribou ~ distHPAR_km*Felev + traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)

AICtab(mfull.1,mfull.2,mfull.3,mfull.4,mfull.5,mfull.0,m0, weights=TRUE) #mfull.1 is the best model by far (weight of 1)


AICtab(mfull.2,mb.3,mh.1,mt.3, weights=TRUE) #mfull.2 is the best model based on AIC and loglikelihood ranking


###--- Top Model
summary(mfull.1)


#do for each park
names(CA.df)
mdata <- CA.df[c(2:5,13,15,21,22,27,29,31)] #smaller dataframe
head(mdata)


CA.NNPR <- mdata %>%
  filter(park=="NNPR")
head(CA.NNPR)

CA.NaNPR <- mdata %>%
  filter(park=="NaNPR")
head(CA.NaNPR)

mfull.NNPR <- glmmTMB(Caribou ~ distHPAR_km*Season +traffic_log + as.factor(Year) + (1|site) + (1|Month), 
                      data=CA.NNPR,family = nbinom2)
mfull.NaNPR <- glmmTMB(Caribou ~ distHPAR_km*Season +traffic_log + as.factor(Year) + (1|site) + (1|Month), 
                       data=CA.NaNPR,family = nbinom2)


summary(mfull.NNPR)
summary(mfull.NaNPR)


#############################################################
###--- Housekeeping - cleaning up files
###--- Saving only camop matrix and camdf covariate dataframe

rm(list= ls()[!(ls() %in% c('CA.df','CA.NNPR','CA.NaNPR','mfull.1','mfull.NNPR','mfull.NaNPR'))])
save.image(file="HPAR_GLMMs.RData") 

###########################################################################
#Family: nbinom2  ( log )
#Formula:          Caribou ~ distHPAR_km * Season + traffic_log + as.factor(Year) +      (1 | site) + (1 | Month)
#Data: CA.df

#mfull.1 <- glmmTMB(Caribou ~ distHPAR_km*Season +traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = nbinom2)
summary(mfull.1)
confint(mfull.1)
exp(confint(mfull.1))
mfull.plot <- confint(mfull.1)
mfull.plot2 <- as.data.frame(mfull.plot[2:13,])
mfull.plot2$Variable <- rownames(mfull.plot2)
mfull.plot2$Park <- "Both"

mfull.plot.NNPR <- confint(mfull.NNPR)
mfull.plot.NNPR2 <- as.data.frame(mfull.plot.NNPR[2:13,])
mfull.plot.NNPR2$Variable <- rownames(mfull.plot.NNPR2)
mfull.plot.NNPR2$Park <- "NNPR"

mfull.plot.NaNPR <- confint(mfull.NaNPR)
mfull.plot.NaNPR2 <- as.data.frame(mfull.plot.NaNPR[2:11,])
mfull.plot.NaNPR2$Variable <- rownames(mfull.plot.NaNPR2)
mfull.plot.NaNPR2$Park <- "NaNPR"
nrow(mfull.plot.NaNPR2)
tail(mfull.plot2)

GLMM.plot <- rbind(mfull.plot2,mfull.plot.NNPR2,mfull.plot.NaNPR2)
tail(GLMM.plot)

# Plot the estimates and 95% CI for each variable 
GLMM.plot
colnames(GLMM.plot)
colnames(GLMM.plot) <- c("CI_2.5","CI_97.5","Estimate","Variable","Park")

GLMM.plot$Variable <- gsub('cond.', '',GLMM.plot$Variable)
levels(as.factor(GLMM.plot$Variable))

GLMM.plot$Variable <- mapvalues(GLMM.plot$Variable,
                                  from = c("as.factor(Year)2014","as.factor(Year)2015","as.factor(Year)2016","as.factor(Year)2017",
                                                                 "distHPAR_km","distHPAR_km:Seasonfall","distHPAR_km:Seasonspring","distHPAR_km:Seasonwinter",
                                                                   "Seasonfall","Seasonspring","Seasonwinter","traffic_log"),
                                  to = c("Year_2014","Year_2015","Year_2016","Year_2017",
                                  "distHPAR_km","distHPAR_km:Season_Fall","distHPAR_km:Season_Spring","distHPAR_km:Season_Winter",
                                  "Season_Fall","Season_Spring","Season_Winter","Traffic (log)"))



p <- ggplot(GLMM.plot,aes(Variable, Estimate)) +
  geom_hline(yintercept = 0, color="darkgrey") +
  geom_point(colour="black", shape=16, size=4)+
  geom_linerange(aes(Variable, ymin = CI_2.5, ymax = CI_97.5))+
  coord_flip()+
  theme(axis.title.y=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14)) +
  facet_grid(.~Park, scales="free")

library(Cairo)
Cairo(file="GLMM_mfull.1.PNG", 
      type="png",
      width=4000, 
      height=2500, 
      pointsize=15,
      bg="white",
      dpi=300)
p
dev.off()


###--- Plotting observed data
levels(CA.df$Season)
CA.df$Season <- factor(CA.df$Season,levels(CA.df$Season)[c(1,2,4,3)])

xyplot(Caribou ~  distHPAR_km| factor(Season),
       data = CA.df, 
       xlab = "Distance to HPAR (km)",
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

xyplot(Caribou ~  traffic_log,
       data = CA.df, 
       xlab = "Traffic volume (log)",
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


boxplot(CA.df$Caribou~as.factor(CA.df$Year), ylab="Monhtly Caribou Detections")
