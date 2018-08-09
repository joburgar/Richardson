#############################################################
# HPAR_JAGS_view.output.R
# created by Joanna Burgar, 17-Mar-2016
# script for viewing JAGS output from folder of output files
#############################################################

# set directory
getwd()
setwd("/Users/JBurgar/Google Drive/Parks Canada - Caribou Contract March 2018/HPAR (1)/Cedar_Output/") # set working directory

######################################
#setwd and load packages
library(coda)
library(xtable)
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
library(Cairo)


######################################
#load files
info <- list.files(getwd(),"mc.*RData")
for(i in 1:length(info)) load(info[[i]]) 

info.mc <- unique(gsub('ET.','',info))
info.mc <- str_sub(info.mc, 1, str_length(info.mc)-6)

#[1] "mc.CA2013F" "mc.CA2014F" "mc.CA2014S" "mc.CA2015F" "mc.CA2015S"

######################################
#create function to load files and write output to table

#function to create output table for JAGS output
get_JAGS_output <- function(filename){
  out <- filename
  s <- summary(window(out, start = 51001))
  gd <- gelman.diag(window(out, start = 51001),multivariate = FALSE)
  output_table <- rbind(as.data.frame(t(s$statistics)),
                        as.data.frame(t(s$quantiles)),
                        as.data.frame(t(gd$psrf)))
  return(output_table)	
}

#function to return model run time (hours)
get_JAGS_ET <- function(filename){
  ET <- parse_number(filename)
  return(ET)
}

######################################


###--- create data.frame for JAGS MRT (same order as loaded)
# load files
info.ET <- list.files(getwd(),"mc.*ET.RData")
info.ET <- str_sub(info.ET, 1, str_length(info.ET)-6)

#[1] "mc.CA2013F.ET" "mc.CA2014F.ET" "mc.CA2014S.ET" "mc.CA2015F.ET" "mc.CA2015S.ET"

info.ET_MRT <- vector('list', length(info.ET))
names(info.ET_MRT) <- paste0('MRT_hrs', seq_along(info.ET))
for(i in seq_along(info.ET_MRT)){
  info.ET.list <- get_JAGS_ET(get(info.ET[i]))
  info.ET_MRT[i] <- info.ET.list
}


MRT <- as.data.frame(unlist(info.ET_MRT))
colnames(MRT) <- c("MRT_hours")
MRT$MRT_days <- MRT$MRT_hours/24
MRT$model <- as.vector(substring(info.ET,4))
MRT$model <- as.factor(str_sub(MRT$model, 1, str_length(MRT$model)-3))
summary(MRT)

write.csv(MRT,file="HPAR_JAGS_MRT.csv")

##################################################################################################
###---format and aggregate model data

JAGS.output <- lapply(info.mc, function(i) get_JAGS_output(get(i)))
JAGS.output <- do.call(rbind, JAGS.output)
glimpse(JAGS.output)

JAGS.output.df <- as.data.frame(JAGS.output)
summary(JAGS.output.df)


JAGS.model <- rep(info.mc, times=1, each=11)
JAGS.output.df$model <- JAGS.model
head(JAGS.output.df,20)

estimates <- c("Mean", "SD", "Naive SE","Time-series SE", "CI_2.5","CI_25", "CI_50", "CI_75", "CI_97.5", "gd.point","gd.upper")
JAGS.output.df$estimates <- as.factor(rep(estimates, nrow(JAGS.output.df)/11))


write.csv(JAGS.output.df, file="HPAR_JAGS_modeloutput.csv")


###--- manipulate data to be in correct format for ggplot graphs
data <- read.csv("HPAR_JAGS_modeloutput.csv") 

glimpse(data)
data <- data[c(-1)]
head(data)

data$Season <- as.factor(substr(data$model,10,10))
data$Year <- as.factor(substr(data$model,6,9))
glimpse(data)
names(data)

data2 <- gather(data, key=parameter, value=value, 1:4)
glimpse(data2)
data2$parameter <- as.factor(data2$parameter)
str(data2)
summary(data2)
table(data2$estimate, data2$parameter)

data.plot <- spread(data2, estimates, value)
head(data.plot)
str(data.plot)

levels(data.plot$parameter)
data.plot$parameter <- factor(data.plot$parameter, levels = c("sigma", "lam0", "psi", "D"))

data2$parameter <- mapvalues(data2$parameter, 
                             from = c("sigma", "lam0", "psi", "D"),
                             to = c("Sigma", "Lam0","Psi","Density"))
data2$parameter <- factor(data2$parameter, levels = c("Density","Lam0","Sigma", "Psi"))

density <- filter(data2, parameter=="Density")
density.plot <- spread(density, estimates, value)
density.plot$model <- as.factor(substr(density.plot$model,6,10))

lam0 <- filter(data2, parameter=="Lam0")
lam0.plot <- spread(lam0, estimates, value)
lam0.plot$model <- as.factor(substr(lam0.plot$model,6,10))

sigma <- filter(data2, parameter=="Sigma")
sigma.plot <- spread(sigma, estimates, value)
sigma.plot$model <- as.factor(substr(sigma.plot$model,6,10))

psi <- filter(data2, parameter=="Psi")
psi.plot <- spread(psi, estimates, value)
psi.plot$model <- as.factor(substr(psi.plot$model,6,10))

data3 <- rbind(density.plot, lam0.plot, sigma.plot, psi.plot)
data3$CV <- data3$SD/data3$Mean
summary(data3$CV)
summary(data3$gd.point) # Gelman-Rubin diagnostic <1.1 suggests model convergence
summary(data3$gd.upper) # Gelman-Rubin diagnostic <1.1 suggests model convergence
write.csv(data3, file="HPAR_JAGS_modeloutput.csv")

##########################################################################
##########################################################################
# Plot the Mean and standard deviation for each model by species:  

data3 <- read.csv("HPAR_JAGS_modeloutput.csv")
head(data3)
data3$Season <- as.factor(data3$Season)
glimpse(data3)
data3$Season <- mapvalues(data3$Season, from = c( "S","F"),to = c("Summer","Fall"))
relevel(data3$Season, ref="Summer")

p = ggplot(data3, aes(model, CI_50))+
  geom_point(colour="white", shape=21, size=10, 
             aes(fill=Season))+
  scale_fill_manual(values=c("orangered4","springgreen4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name="Mode ± BCI") +
  geom_linerange(aes(model, ymin = CI_2.5, ymax = CI_97.5)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~parameter, scales="free_y")

###--- create graph for only Density and Lam0 (probability of encounter/encounter rate)
p2 = ggplot(data3[grep("Density|Lam0", data3$parameter),], aes(model, CI_50))+
  geom_point(colour="white", shape=21, size=10, 
             aes(fill=Season))+
  scale_fill_manual(values=c("orangered4","springgreen4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name="Mode ± BCI") +
  geom_linerange(aes(model, ymin = CI_2.5, ymax = CI_97.5)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.y = element_text(size=14))+
  theme(axis.text.x = element_text(size=14))+
  facet_wrap(~parameter, scales="free_y")

Cairo(file="HPAR_D&Lam0_output.PNG", 
      type="png",
      width=3000, 
      height=2000, 
      pointsize=15,
      bg="white",
      dpi=300)
p2
dev.off()

