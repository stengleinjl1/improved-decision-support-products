# Leveraging modern quantitative tools to improve wildlife decision support products for natural resources agencies
# Author: Jennifer Stenglein

# Vignette 2: Population model inputs from nonparametric spatial regression
# Vignette 3: Population model outputs from a Bayesian model

# Code to create and run a Bayesian model that uses all available data and nonparametric spatial regression 
# to estimate yearling buck percent, yearling doe percent, and fawn-to-doe ratio for white-tailed deer by deer management unit
# in the same Bayesian model, a sex-age-kill population estimate is constructed as a derived quantity.
# Example is for Wisconsin, USA in 2022

library(tidyverse)
library(jagsUI)

# read in raw data
dat <- read_csv("data/dmudata2022.csv")
dat <- dat %>%
  mutate(forest.binary = ifelse(zone == "Forest",1,0))

# read in knots
# 50 knots were found with cover.design package and read in here for simplicity
knots <- read_csv("data/knots.csv")

# scale coordinates 
mean_x <- mean(dat$X)
sd_x <- sd(dat$X)
mean_y <- mean(dat$Y)
sd_y <- sd(dat$Y)

knots <- knots %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y) %>%
  dplyr::select(X.scale,Y.scale)

dat <- dat %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y)

# get matrix ready for spatial smoothing
knots.dist <- dist(knots,"euclidean",diag=T,upper=T)
omega_all = knots.dist^2*log(knots.dist) # basis
svd.omega_all <- svd(omega_all)
sqrt.omega_all <- t(svd.omega_all$v %*%
                      (t(svd.omega_all$u)*sqrt(svd.omega_all$d)))

# further prepare data for model input
# fdr data prep - keep only locations that have data to contribute
dat.fdr <- dat %>%
  dplyr::select(DMU,forest.binary,X.scale,Y.scale,fawns.sdo,does.sdo) %>%
  mutate(total.sdo = fawns.sdo + does.sdo) %>%
  filter(total.sdo>0) %>%
  arrange(DMU)
n.fdr = dim(dat.fdr)[1] # number of locations with data

# ybp data prep - keep only locations that have data to contribute
dat.ybp <- dat %>%
  dplyr::select(DMU,forest.binary,X.scale,Y.scale,bucks.aged.yearling,bucks.aged.adult) %>%
  filter(bucks.aged.adult>0) %>%
  arrange(DMU)
n.ybp = dim(dat.ybp)[1] # number of locations with data

# ydp data prep - keep only locations that have data to contribute
dat.ydp <- dat %>%
  dplyr::select(DMU,forest.binary,X.scale,Y.scale,does.aged.yearling,does.aged.adult) %>%
  filter(does.aged.adult>0) %>%
  arrange(DMU)
n.ydp = dim(dat.ydp)[1] # number of locations with data

# prepare prediction dataset to include all DMUs
dat.pred = dat %>%
  dplyr::select(DMU,forest.binary,X.scale,Y.scale,harv.total,harv.antlered) %>%
  na.omit() %>%
  distinct() %>%
  arrange(DMU)
N = dim(dat.pred)[1] # number of predicted locations (DMU locations)

# prediction dataset for spatial smoothing
cov.dist_pred = fields::rdist(x1=cbind(dat.pred$X.scale,dat.pred$Y.scale),x2=knots)
Z_K_pred = cov.dist_pred^2*log(cov.dist_pred)
Z.pred <- t(solve(sqrt.omega_all,t(Z_K_pred)))

# standardize for better performance
meanZ = mean(Z.pred)
sdZ = sd(Z.pred)
Z.pred <- (Z.pred - meanZ)/sdZ # this should still be mean(Z) so that all and pred are on the same scale

# now for each spline
cov.dist_all.fdr = fields::rdist(x1=cbind(dat.fdr$X.scale,dat.fdr$Y.scale),x2=knots)
Z_K.fdr = cov.dist_all.fdr^2*log(cov.dist_all.fdr) # basis
Z.fdr <- t(solve(sqrt.omega_all,t(Z_K.fdr)))
Z.fdr <- (Z.fdr - meanZ)/sdZ 

cov.dist_all.ybp = fields::rdist(x1=cbind(dat.ybp$X.scale,dat.ybp$Y.scale),x2=knots)
Z_K.ybp = cov.dist_all.ybp^2*log(cov.dist_all.ybp) # basis
Z.ybp <- t(solve(sqrt.omega_all,t(Z_K.ybp)))
Z.ybp <- (Z.ybp - meanZ)/sdZ 

cov.dist_all.ydp = fields::rdist(x1=cbind(dat.ydp$X.scale,dat.ydp$Y.scale),x2=knots)
Z_K.ydp = cov.dist_all.ydp^2*log(cov.dist_all.ydp) # basis
Z.ydp <- t(solve(sqrt.omega_all,t(Z_K.ydp)))
Z.ydp <- (Z.ydp - meanZ)/sdZ 

# bundle data for model run
jags.data = list(
  fdr.fawn.N = dat.fdr$fawns.sdo,
  fdr.doe.N = dat.fdr$does.sdo,
  fdr.fawndoe.N = dat.fdr$total.sdo,
  fdr.forest = dat.fdr$forest.binary,
  Z.fdr = Z.fdr,
  fdr.ndmus = n.fdr,
  
  ybp.bucks.yearlings.N = dat.ybp$bucks.aged.yearling,
  ybp.bucks.adults.N = dat.ybp$bucks.aged.adult,
  ybp.forest = dat.ybp$forest.binary,
  Z.ybp = Z.ybp,
  ybp.ndmus = n.ybp,
  
  ydp.does.yearlings.N = dat.ydp$does.aged.yearling,
  ydp.does.adults.N = dat.ydp$does.aged.adult,
  ydp.forest = dat.ydp$forest.binary,
  Z.ydp = Z.ydp,
  ydp.ndmus = n.ydp,
  
  forest.pred = dat.pred$forest.binary,
  Z.pred = Z.pred,
  buckharv = dat.pred$harv.antlered,
  totalharv = dat.pred$harv.total,
  
  nknots = 50,
  Ndmus = N
)

# initial values
inits <- function() list(
  beta.fdr.fawn = runif(2,-1,1),
  beta.fdr.doe = runif(2,-1,1),
  beta.ydp = runif(2,-1,1),
  beta.ybp = runif(2,-1,1)
)

# parameters to track
params <- c("beta.fdr.fawn",
           "beta.fdr.doe",
           "beta.ydp",
           "beta.ybp",
           "sigma.fdr.fawn.spline",
           "sigma.fdr.doe.spline",
           "sigma.ydp.spline",
           "sigma.ybp.spline",
           "fdr",
           "ydp",
           "ybp",
           "posthunt")

# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 3

# model
sink("model.txt")
cat("
model {
  
# PRIORS 
  
  # intercept and binary forest priors
  for (p in 1:2) {
    beta.fdr.fawn[p] ~ dunif(-2,2)
    beta.fdr.doe[p] ~ dunif(-2,2)
    beta.ydp[p] ~ dunif(-2,2)
    beta.ybp[p] ~ dunif(-2,2)
  }
 
  # spline random effect priors
  sigma.fdr.fawn.spline~dunif(0,100)
  tau.fdr.fawn.spline <- 1/(sigma.fdr.fawn.spline*sigma.fdr.fawn.spline)
  sigma.fdr.doe.spline~dunif(0,100)
  tau.fdr.doe.spline <- 1/(sigma.fdr.doe.spline*sigma.fdr.doe.spline)
  sigma.ydp.spline~dunif(0,100)
  tau.ydp.spline <- 1/(sigma.ydp.spline*sigma.ydp.spline)
  sigma.ybp.spline~dunif(0,100)
  tau.ybp.spline <- 1/(sigma.ybp.spline*sigma.ybp.spline)

  for (k in 1:nknots) {
    alpha.fdr.fawn.spline[k] ~ dnorm(0,tau.fdr.fawn.spline)
    alpha.fdr.doe.spline[k] ~ dnorm(0,tau.fdr.doe.spline) 
    alpha.ydp.spline[k] ~ dnorm(0,tau.ydp.spline)
    alpha.ybp.spline[k] ~ dnorm(0,tau.ybp.spline)
  }
  
  # prior for buck recovery rate
  brr ~ dunif(0.65,0.75)

# LIKELIHOOD

  for (i in 1:fdr.ndmus){
   
    fdr.fawn.N[i] ~ dbin(fdr.fawn.p[i],fdr.fawndoe.N[i])
    fdr.doe.N[i] ~ dbin(fdr.doe.p[i],fdr.fawndoe.N[i])
    
    logit(fdr.fawn.p[i]) <- beta.fdr.fawn[1] + beta.fdr.fawn[2]*fdr.forest[i] + fdr.fawn.spline[i]
    logit(fdr.doe.p[i]) <- beta.fdr.doe[1] + beta.fdr.doe[2]*fdr.forest[i] + fdr.doe.spline[i]
    
    fdr.fawn.spline[i] <- inprod(alpha.fdr.fawn.spline[1:nknots],Z.fdr[i,1:nknots])
    fdr.doe.spline[i] <- inprod(alpha.fdr.doe.spline[1:nknots],Z.fdr[i,1:nknots])
  }

  for (i in 1:ydp.ndmus) {
  
    ydp.does.yearlings.N[i] ~ dbin(ydp.p[i],ydp.does.adults.N[i])
    
    logit(ydp.p[i]) <- beta.ydp[1] + beta.ydp[2]*ydp.forest[i] + ydp.spline[i]
    
    ydp.spline[i] <- inprod(alpha.ydp.spline[1:nknots],Z.ydp[i,1:nknots])
  
  }

  for (i in 1:ybp.ndmus) {
  
    ybp.bucks.yearlings.N[i] ~ dbin(ybp.p[i],ybp.bucks.adults.N[i])
    
    logit(ybp.p[i]) <- beta.ybp[1] + beta.ybp[2]*ybp.forest[i] + ybp.spline[i]
    
    ybp.spline[i] <- inprod(alpha.ybp.spline[1:nknots],Z.ybp[i,1:nknots])
  
  }

# PREDICTION

  for (i in 1:Ndmus) {
  
    logit(fdr.fawn.p.pred[i]) <- beta.fdr.fawn[1] + beta.fdr.fawn[2]*forest.pred[i] + fdr.fawn.spline.pred[i]
    logit(fdr.doe.p.pred[i]) <- beta.fdr.doe[1] + beta.fdr.doe[2]*forest.pred[i] + fdr.doe.spline.pred[i]
    
    logit(ydp[i]) <- beta.ydp[1] + beta.ydp[2]*forest.pred[i] + ydp.spline.pred[i]
    
    logit(ybp[i]) <- beta.ybp[1] + beta.ybp[2]*forest.pred[i] + ybp.spline.pred[i]
    
    fdr.fawn.spline.pred[i] <- inprod(alpha.fdr.fawn.spline[1:nknots],Z.pred[i,1:nknots])
    fdr.doe.spline.pred[i] <- inprod(alpha.fdr.doe.spline[1:nknots],Z.pred[i,1:nknots])
    
    ydp.spline.pred[i] <- inprod(alpha.ydp.spline[1:nknots],Z.pred[i,1:nknots])
    
    ybp.spline.pred[i] <- inprod(alpha.ybp.spline[1:nknots],Z.pred[i,1:nknots])
    
    # DERIVED QUANTITIES
    fdr[i] <- fdr.fawn.p.pred[i]/fdr.doe.p.pred[i]
    
    # adult does to adult bucks ratio (not including the difference in sex ratio - ydp*1.1)
    dbr[i] <- ybp[i] / ydp[i]
    
    bucks[i] <- buckharv[i] /  (ybp[i]*brr)
    does[i] <- bucks[i] * dbr[i]
    fawns[i] <- does[i] * fdr[i]
    
    prehunt[i] <- bucks[i]+does[i]+fawns[i]
    posthunt[i] <- prehunt[i] - (totalharv[i]*1.15)
    
  }
}
", fill = T)
sink()

model <- read_file("model.txt")

# Run model
out = jags(data = jags.data, 
           inits = inits, 
           parameters.to.save = params, 
           model = "model.txt", 
           n.chains = nc, 
           n.iter = ni, 
           n.burn = nb, 
           n.thin = nt, 
           parallel = T)
