##############################################################
# HPAR_PowerAnalysis.R
# Created by Joanna Burgar, 015-Mar-2018
# Running power analysis for GLMMs
#############################################################

###--- Load packages
library(simr) # power analysis package
#library(lme4) # glmm package associated with simr (should be loaded if loaded simr with dependencies)

###--- Load data - 
setwd("C:/Users/JBurgar/Documents/R/Analysis/HPAR/") # set working directory
load("HPAR_GLMMs.RData")

###--- run glmer models - non-convergence with interaction but as primarily interested in distHPAR_km, remove interaction for power analysis
###--- lme4 does not have capabilities for negative binomial so using poisson distribution instead

p.NNPR <- glmer(Caribou ~ distHPAR_km + Season + traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.NNPR,family = poisson)
p.NaNPR <- glmer(Caribou ~ distHPAR_km + Season + traffic_log + (1|site) + (1|Month), data=CA.NaNPR,family = poisson) # removed Year as would not converge
p.full <- glmer(Caribou ~ distHPAR_km + Season + traffic_log + as.factor(Year) + (1|site) + (1|Month), data=CA.df,family = poisson)

###################################################################################
#run power first for full model (both parks)
summary(p.full) # estimated effect size for distHPAR_km is -0.07788323, is NOT sig @ 0.05 level using default t-test
fixef(p.full)["distHPAR_km"] 
#distHPAR_km 
#-0.07788323 

#consider the power required to detect a slope of -0.5
#fixef(p.full)["distHPAR_km"]  <- -0.5
#ps.p.full1 <- powerSim(p.full)

###--- provides the percentage for the power to reject the null hypothesis of zero trend in x 
#Power for predictor 'distHPAR_km', (95% confidence interval):
#  100.0% (99.63, 100.0) # 80% power is considered adequate - 100% is more than enough power
#Test: z-test
#Effect size for distHPAR_km is -1.0
#Based on 1000 simulations, (59 warnings, 0 errors) - 59/100 simulations where data was a poor fit for the model
#alpha = 0.05, nrow = 3055
#Time elapsed: 1 h 4 m 10 s


###################################################################################
#run power for NaNPR & NNPR
summary(p.NaNPR) # estimated effect size for distHPAR_km is 0.12723, is NOT sig @ 0.05 level using default t-test
fixef(p.NaNPR)["distHPAR_km"] 
#distHPAR_km 
#0.127226

#consider the power required to detect a slope of -0.5
fixef(p.NaNPR)["distHPAR_km"]  <- -0.5
ps.NaNPR <- powerSim(p.NaNPR)

###--- provides the percentage for the power to reject the null hypothesis of zero trend in x 


summary(p.NNPR) # estimated effect size for distHPAR_km is -0.17448, is NOT sig @ 0.05 level using default t-test
fixef(p.NNPR)["distHPAR_km"] 
#distHPAR_km 
#-0.1744758

#consider the power required to detect a slope of -1.0
fixef(p.NNPR)["distHPAR_km"]  <- -0.5
ps.NNPR <- powerSim(p.NNPR)


###################################################################################
###--- adding more sites

length(unique(CA.NaNPR$site)) # NaNPR - 17 sites
length(unique(CA.NNPR$site)) # NNPR - 44 sites
# add another 23 & 21 sites for a total of 40 sites for NaNPR and 65 sites for NNPR

#calculate the effect of increasing the number of sites to 40 for NaNPR
p.NaNPR.site <- extend(p.NaNPR, along="site", n=40)
pc.NaNPR.site <- powerCurve(p.NaNPR.site, along="site")

png('pc.NaNPR.site.png')
plot(pc.NaNPR.site)
dev.off()

#calculate the effect of increasing the number of sites to 65 for NNPR
p.NNPR.site <- extend(p.NNPR, along="site", n=65)
pc.NNPR.site <- powerCurve(p.NNPR.site, along="site")

png('pc.NNPR.site.png')
plot(pc.NNPR.site)
dev.off()



###################################################################################
###--- adding more months

table(CA.NaNPR$Yr_Month, CA.NaNPR$Caribou) # NanPR first caribou detection in July 2015, last caribou detection in June 2016
head(table(CA.NNPR$Yr_Month, CA.NNPR$Caribou))   # NNPR first caribou detection in Aug 2013
tail(table(CA.NNPR$Yr_Month, CA.NNPR$Caribou))   # NNPR last caribou detection in June 2016

length(unique(CA.NaNPR$Yr_Month)) #47 months (JUly 2013 - July 2017)
length(unique(CA.NNPR$Yr_Month)) #44 months
# extend for another 2 years (24 months) for a total of 72 months for both parks

# calculate the effect of increasing the number of study months to 72 (6 years)
# NaNPR
p.NaNPR.month <- extend(p.NaNPR, along="Yr_Month", n=72)
pc.NaNPR.month <- powerCurve(p.NaNPR.month , along="Yr_Month")

png('pc.NaNPR.month.png')
plot(pc.NaNPR.month)
dev.off()

# NNPR
p.NNPR.month <- extend(p.NNPR, along="Yr_Month", n=72)
pc.NNPR.month <- powerCurve(p.NNPR.month , along="Yr_Month")

png('pc.NNPR.month.png')
plot(pc.NNPR.month)
dev.off()

save.image(file="HPAR_PowAnalysis.RData") 
#load("HPAR_PowAnalysis.RData")
ls()
