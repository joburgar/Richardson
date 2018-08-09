#############################################################
# 8.RICH_SC_cloud.R
# SCR book Ch 18 code;
# created by Joanna Burgar, 09-Aug-2018 for use on compute canada
# script for estimating density of unmarked caribou population
# using JAGS with vaguely informed priors
#############################################################

###--- set directory and load files

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

# specify how much to buffer sites by (in 10km units, according to coord.scale)
coord.scale <- 10000
buffer <- 10 # 10 km unit buffer

traplocs <- as.matrix(traplocs)           # 30 camera locations in this matrix
X <- traplocs/coord.scale                     # 30 trap coords in 10 km units
dim(X)

###--- format camop to operability matrix for use in models
# remove sites 5,12,13,24,26-28,35,41,42,45,49-65
oper <- camop
oper[is.na(oper)] <- 0

###--- create xlims and ylims of scaled coordinates
## SnowFree
X.scaled <- X[,1]-min(X[,1])

Xl.scaled <- min(X.scaled - buffer)
Xu.scaled <- max(X.scaled + buffer)

Y.scaled <- X[,2]-min(X[,2])
Yl.scaled <- min(Y.scaled - buffer)
Yu.scaled <- max(Y.scaled + buffer)

xlims.scaled <- c(Xl.scaled,Xu.scaled); ylims.scaled <- c(Yl.scaled,Yu.scaled)  

areakm2.scaled <- (Xu.scaled - Xl.scaled)*(Yu.scaled - Yl.scaled)*100
# [1] 1838 km2

X2 <- as.matrix(cbind(X.scaled,Y.scaled))
dim(X2) # scaled traploc matrix in 10 km units

###--- remaining model parameters
500/1838 # density of 0.27 caribou/km2

#############################################################
###--- RUNNING 3 CHAINS IN PARALLEL USING JAGS IMPLEMENTATION

dat <- list(n=n, X=X2, M=500, J=nrow(n), K=ncol(n), xlim=xlims.scaled, ylim=ylims.scaled, area=areakm2.scaled)

# specify initial values
init <-  function() {  list(sigma=rnorm(1,10), lam0=runif(1) , z=rep(1,dat$M)) }

# specify parameters to monitor
pars <- c("sigma","lam0","psi","D")


# specify model
cat("
    model {
    sigma ~ dgamma(15,50) # using home range area (38-620 km2) to form vaguely informed priors
    lam0 ~ dunif(0,10)
    psi ~ dbeta(1,1)
    for(i in 1:M) {
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1], xlim[2])
    s[i,2] ~ dunif(ylim[1], ylim[2])
    for(j in 1:J) {
    distsq[i,j] <- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
    lam[i,j] <- lam0 * exp(-distsq[i,j] / (2*sigma^2)) * z[i]
    }
    }
    for(j in 1:J) {
    bigLambda[j] <- sum(lam[,j])
    for(k in 1:K) {
    n[j,k] ~ dpois(bigLambda[j])
    }
    }
    N <- sum(z[])
    D <- N/area
    }
    ",file="CAJ.txt")

# run model in jags (and calculate model run time)
#Snow
library(rjags)
library(coda)
library(parallel)

(start.time <- Sys.time())
cl3 <- makeCluster(3)
clusterExport(cl3, c("dat","init","pars"))
CA.JAGS <- clusterEvalQ(cl3, {
  library(rjags)
  jm2 <- jags.model("CAJ.txt", data=dat, inits=init,
                    n.chains=1, n.adapt=1000)
  jc2 <- coda.samples(jm2, pars, n.iter=100000)
  return(as.mcmc(jc2))
})
mc.CA2013F <- mcmc.list(CA.JAGS)

(end.time <- Sys.time())
mc.CA2013F.ET <- difftime(end.time, start.time, units='hours')

save("mc.CA2013F.ET",file="mc.CA2013F.ET.RData")
save("mc.CA2013F",file="mc.CA2013F.RData")
