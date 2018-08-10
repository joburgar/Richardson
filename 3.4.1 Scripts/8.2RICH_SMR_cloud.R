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

traplocs <- camdf[,2:3]           # 65 camera locations in this matrix
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
    sigma ~ dgamma(24,8) # weakly informative prior - 38-619 km2
    lam0 ~ dunif(0,10)
    psi ~ dbeta(1,1)
    
    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    
    for(j in 1:J) {
    distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
    lam[i, j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i] 
    }
    } # End of 1:M
    for(j in 1:J) {
    for(k in 1:K) {
    bigLambda[j, k] <- sum(lam[1:M, j]) * oper[j, k]
    n[j, k] ~ dpois(bigLambda[j, k])
    
    }
    }
    N <- sum(z[])
    D <- N/area
    }
    ",file="CA.txt")

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

