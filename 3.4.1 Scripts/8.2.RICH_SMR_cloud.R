#!/bin/env Rscript

#############################################################
# 8.RICH_SC_cloud.R
# SCR book Ch 18 code;
# created by Joanna Burgar, 09-Aug-2018 for use on compute canada
# script for estimating density of unmarked caribou population
# using JAGS with vaguely informed priors
#############################################################

###--- set directory and load files
setwd("./") # for compute canada

setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.2 Input") # set working directory

###--- Load data

load("camop.RData") # load camera effort data
load("SpMatrix.RData")# load species detection data

traplocs <- camdf[,2:3]           # 30 camera locations in this matrix
colnames (traplocs) <- c("x","y")
str(traplocs)

###--- Caribou detectopms
dim(CA.matrix) # 30 x 151

# 2017- snow, maximising collared detections
sum(CA.matrix[,4:68],na.rm=T)  #21 detections 14 Nov - 13 Jan
CA.2017SN <- CA.matrix[,4:68]

n <- as.matrix(CA.2017SN)
n[is.na(n)] <- 0 # changes all NA, as camop will deal with inactive cameras
sum(n)

###--- create encounter history data frame for the 2 caribou detected on 3 days
DayLookup[which(DayLookup$Datep=="2017-11-15"),] # 5
DayLookup[which(DayLookup$Datep=="2017-11-23"),] # 13
DayLookup[which(DayLookup$Datep=="2017-12-17"),] # 37

edf <- matrix(NA,3,4)
edf <- as.data.frame(edf)
colnames(edf) <- c("Sess","Indiv","Occ","Trap")
edf$Sess <- 1
edf$Indiv <- c("2607","2607","2228")
edf$Occ <- c(DayLookup[which(DayLookup$Datep=="2017-11-15"),1],
             DayLookup[which(DayLookup$Datep=="2017-11-23"),1],
             DayLookup[which(DayLookup$Datep=="2017-12-17"),1])
edf$Trap <- c(19,24,20)

str(edf)

###---- camera station coordinates
summary(traplocs)

# specify how much to buffer sites by (in 1 km units, according to coord.scale)
coord.scale <- 1000
buffer <- 10 # 10 km unit buffer

traplocs <- as.matrix(traplocs)           # 30 camera locations in this matrix
X <- traplocs/coord.scale                     # 30 trap coords in 1 km units
dim(X)

###--- create xlims and ylims of scaled coordinates
X.scaled <- X[,1]-min(X[,1])

Xl.scaled <- min(X.scaled - buffer)
Xu.scaled <- max(X.scaled + buffer)

Y.scaled <- X[,2]-min(X[,2])
Yl.scaled <- min(Y.scaled - buffer)
Yu.scaled <- max(Y.scaled + buffer)

xlims.scaled <- c(Xl.scaled,Xu.scaled); ylims.scaled <- c(Yl.scaled,Yu.scaled)  

areakm2.scaled <- (Xu.scaled - Xl.scaled)*(Yu.scaled - Yl.scaled)
# [1] 1736 km2

X2 <- as.matrix(cbind(X.scaled,Y.scaled))
dim(X2) # scaled traploc matrix in 1 km units

###--- remaining model parameters
500/1736 # density of 0.29 caribou/km2
200/1736 # density of 0.12 caribou/km2

M <- 200

###--- format camop to operability matrix for use in models
oper <- camop
oper[is.na(oper)] <- 0
oper <- oper[,4:68]
dim(oper) # 30 traps by 65 days

###--- trap operability calculations to determine the number of active traps and days out of potential
#sum(oper) #4231 possible camera days over entire time; 1844 during suvery period
#sum(rowSums(oper!= 0))
#Roper <- rowSums(oper)
#Roper[order(Roper)]
#sum(Roper!= 0) # 30 sites operational between 16-151 days; 15-65 for survey period
#Coper <- colSums(oper)
#Coper[order(Coper)] # 10 - 30 sites active each day; between 26-30 for survey period


#############################################################
###--- RUNNING 3 CHAINS IN PARALLEL USING JAGS IMPLEMENTATION
dat <- list(n=n, X=X2, M=M, J=nrow(n), K=ncol(n), xlim=xlims.scaled, ylim=ylims.scaled, area=areakm2.scaled, oper=oper)

# specify initial values
init <-  function() {  list(sigma=rnorm(1,10), lam0=runif(1) , z=rep(1,dat$M)) }

# specify parameters to monitor
pars <- c("sigma","lam0","psi","D","N")


############################################################

###--- Model = Weakly Informative prior

# specify model
cat("
    model {
    
    lam0.mark ~ dunif(0,10) # Baseline encounter rate for marking occasions - same as SC models
    lam0.resight ~ dunif(0,2) # Resighting occasions - same as Jesse's code
    sigma ~ dgamma(24,8) # weakly informative prior - 38-619 km2 (same as for SC models)
    #sigma ~ dunif(0,2) # from Jesse's code
    sigma2 <- sigma^2
    psi ~ dbeta(1,1) # same as SC code
    
    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])

    for(j in 1:J.mark) {   # Marking observations
    d.mark2[i,j] <- (s[i,1]-x.mark[j,1])^2 + (s[i,2]-x.mark[j,2])^2
    lambda.mark[i,j] <- lam0.mark*exp(-d.mark2[i,j]/(2*sigma2))*z[i]
    y.mark[i,j] ~ dbinom(lambda.mark[i,j], K.mark) # Encounter histories from capture process for all individuals including augmented population.
    }

    # Resight - Calculate Encounter Rates for all individuals and detector
    for(j in 1:J.resight) {
    d.resight2[i,j] <- (s[i,1]-x.resight[j,1])^2 + (s[i,2]-x.resight[j,2])^2
    lambda.resight[i,j] <- lam0.resight * exp(-d.resight2[i,j]/(2*sigma2))*z[i]
    }
    }
    
    for (i in 1:n.marked){
    for (j in 1:J.resight){
    for (k in 1:K.resight){
    y.marked.resight[i,j,k] ~ dpois(lambda.resight[i,j])
    }
    }
    }
    for (i in 1:M){
    for (j in 1:J.resight){
    for (k in 1:K.resight){
    y.l[i,j,k] ~ dpois(lambda.resight[i,j]*(1 - marked[i]))
    }
    }
    }
    # Resight likelihoods of unmarked animals
    for (j in 1:J.resight){
    for (k in 1:K.resight){
    n.unmarked[j,k] ~ dsum(y.l[1,j,k], y.l[2,j,k], y.l[3,j,k], y.l[4,j,k], y.l[5,j,k], y.l[6,j,k], y.l[7,j,k], y.l[8,j,k], y.l[9,j,k], y.l[10,j,k], y.l[11,j,k], y.l[12,j,k], y.l[13,j,k], y.l[14,j,k], y.l[15,j,k], y.l[16,j,k], y.l[17,j,k], y.l[18,j,k], y.l[19,j,k], y.l[20,j,k], y.l[21,j,k], y.l[22,j,k], y.l[23,j,k], y.l[24,j,k], y.l[25,j,k], y.l[26,j,k], y.l[27,j,k], y.l[28,j,k], y.l[29,j,k], y.l[30,j,k], y.l[31,j,k], y.l[32,j,k], y.l[33,j,k], y.l[34,j,k], y.l[35,j,k], y.l[36,j,k], y.l[37,j,k], y.l[38,j,k], y.l[39,j,k], y.l[40,j,k], y.l[41,j,k], y.l[42,j,k], y.l[43,j,k], y.l[44,j,k], y.l[45,j,k], y.l[46,j,k], y.l[47,j,k], y.l[48,j,k], y.l[49,j,k], y.l[50,j,k], y.l[51,j,k], y.l[52,j,k], y.l[53,j,k], y.l[54,j,k], y.l[55,j,k], y.l[56,j,k], y.l[57,j,k], y.l[58,j,k], y.l[59,j,k], y.l[60,j,k], y.l[61,j,k], y.l[62,j,k], y.l[63,j,k], y.l[64,j,k], y.l[65,j,k], y.l[66,j,k], y.l[67,j,k], y.l[68,j,k], y.l[69,j,k], y.l[70,j,k], y.l[71,j,k], y.l[72,j,k], y.l[73,j,k], y.l[74,j,k], y.l[75,j,k], y.l[76,j,k], y.l[77,j,k], y.l[78,j,k], y.l[79,j,k], y.l[80,j,k], y.l[81,j,k], y.l[82,j,k], y.l[83,j,k], y.l[84,j,k], y.l[85,j,k], y.l[86,j,k], y.l[87,j,k], y.l[88,j,k], y.l[89,j,k], y.l[90,j,k], y.l[91,j,k], y.l[92,j,k], y.l[93,j,k], y.l[94,j,k], y.l[95,j,k], y.l[96,j,k], y.l[97,j,k], y.l[98,j,k], y.l[99,j,k], y.l[100,j,k], y.l[101,j,k], y.l[102,j,k], y.l[103,j,k], y.l[104,j,k], y.l[105,j,k], y.l[106,j,k], y.l[107,j,k], y.l[108,j,k], y.l[109,j,k], y.l[110,j,k], y.l[111,j,k], y.l[112,j,k], y.l[113,j,k], y.l[114,j,k], y.l[115,j,k], y.l[116,j,k], y.l[117,j,k], y.l[118,j,k], y.l[119,j,k], y.l[120,j,k], y.l[121,j,k], y.l[122,j,k], y.l[123,j,k], y.l[124,j,k], y.l[125,j,k], y.l[126,j,k], y.l[127,j,k], y.l[128,j,k], y.l[129,j,k], y.l[130,j,k], y.l[131,j,k], y.l[132,j,k], y.l[133,j,k], y.l[134,j,k], y.l[135,j,k], y.l[136,j,k], y.l[137,j,k], y.l[138,j,k], y.l[139,j,k], y.l[140,j,k], y.l[141,j,k], y.l[142,j,k], y.l[143,j,k], y.l[144,j,k], y.l[145,j,k], y.l[146,j,k], y.l[147,j,k], y.l[148,j,k], y.l[149,j,k], y.l[150,j,k],
                            y.l[151,j,k],y.l[152,j,k],y.l[153,j,k],y.l[154,j,k],y.l[155,j,k],)
    }
    } 
    
    # Telemetry data for n.collared animals
    for (i in 1:n.collar){
    for (j in 1:n.locs){
    telemetry.array[i, j, 1] ~ dnorm(s[i, 1], 1 / sigma2)
    telemetry.array[i, j, 2] ~ dnorm(s[i, 2], 1 / sigma2)
    }
    }
    
    N <- sum(z[1:M])
    }
    ",file="CA_gSMR.txt")

# run model in jags (and calculate model run time)
library(rjags)
library(coda)
library(parallel)

(start.time <- Sys.time())
cl3 <- makeCluster(3)
clusterExport(cl3, c("dat","init","pars"))
CA.JAGS <- clusterEvalQ(cl3, {
  library(rjags)
  jm2 <- jags.model("CA.txt", data=dat, inits=init, n.chains=1, n.adapt=1000)
  jc2 <- coda.samples(jm2, pars, n.iter=80000)
  return(as.mcmc(jc2))
})
mc.CA2_M200 <- mcmc.list(CA.JAGS)

(end.time <- Sys.time())
mc.CA2_M200.ET <- difftime(end.time, start.time, units='hours')

setwd("C:/Users/JBurgar/Google Drive/Richardson Wildfire Project/3. Data/3.4 Data Analysis/3.4.3 Output") # set ouptput directory
save("mc.CA2_M200",file="mc.CA2_M200.RData")
save("mc.CA2_M200.ET",file="mc.CA2_M200.ET.RData")
stopCluster(cl3)

