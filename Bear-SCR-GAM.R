# ============================================================
# Spatial Capture–Recapture (SCR) Model for Cantabrian Brown Bear
# ============================================================
# Author: José Jiménez
# Institution: Instituto de Investigación en Recursos Cinegéticos (IREC, CSIC-UCLM-JCCM)
# Contact: Jose.Jimenez@csic.es
# ORCID: https://orcid.org/0000-0003-0607-6973
#
# Description:
# This script implements a Bayesian Spatial Capture–Recapture (SCR) model
# to estimate the population size and spatial density of the Cantabrian brown bear (Ursus arctos)
# across its entire range in northern Spain. The model integrates:
#   - Generalized Additive Models (GAMs) for flexible spatial density estimation
#   - Finite mixture structures to account for individual-level movement heterogeneity
#   - Multiple non-invasive genetic sampling methods (scat and hair)
#
# Objectives:
#   1. Estimate total population size and sex ratio
#   2. Identify spatial density patterns and population cores
#   3. Model individual movement behavior using latent mixture groups
#   4. Inform conservation planning by capturing ecological realism
#
# Key Features:
#   - Density modeled as an inhomogeneous Poisson process via spatial splines
#   - Detection modeled using Poisson encounter rates 
#   - Movement modeled with sex-specific finite mixtures
#   - Integrated data from transects and hair traps across four regions
#
# Software:
#   - R (packages: secr, jagam, nimble, coda)
#   - BUGS language for Bayesian inference
#
# Data:
#   - Genetic samples from 16,700 km² study area
#   - Sampling effort and detection histories aligned across methods
#
# Output:
#   - Posterior estimates of population size, density surfaces, movement parameters
#   - Maps of realized density and individual activity centers
#   - Model diagnostics and convergence statistics
#
# Reference:
#   Jiménez et al. (2025). Flexible spatial modelling improves population 
#   estimates for elusive carnivores in fragmented landscapes.
# ============================================================

setwd('C:/...')

library(nimble)
nimble:::setNimbleOption('allowDynamicIndexing', TRUE)
library(basicMCMCplots)
library(coda)
library(lattice)
library(raster)
library(secr)
library(jagsUI)
library(makeJAGSmask)
library(terra)
library(scrbook)
library(mcmcOutput)
library(MCMCvis)
library(tidyverse)
library(ggplot2)

source('SCR_bear_functions.R')

# Scat sampling
bearS.ch <- read.capthist("./data/captScat.txt", "./data/traps.txt", detector='count', noccasions=1)
bearS<-aperm(bearS.ch,c(1,3,2))
yS<-apply(bearS,c(1,2),sum)
(nind<-dim(yS)[1])
yS[,1]<-0
sum(yS)
dim(yS[apply(yS,1,sum)>0,])

# Hair sampling
bearH.ch <- read.capthist("./data/captHair.txt", "./data/traps.txt", detector='count', noccasions=1)
bearH<-aperm(bearH.ch,c(1,3,2))
yH<-apply(bearH,c(1,2),sum)
(nind<-dim(yH)[1])
yH[,1]<-0
sum(yH)
dim(yH[apply(yH,1,sum)>0,])

# Hair-trap sampling Castilla y León
bearHTCyL.ch <- read.capthist("./data/captHairTrapsCyL.txt", "./data/trapPeloCyL.txt", detector='count', noccasions=4)
bearHTCyL<-aperm(bearHTCyL.ch,c(1,3,2))
yHTCyL<-apply(bearHTCyL,c(1,2),sum)
(nind<-dim(yHTCyL)[1])
yHTCyL[,100]<-0
sum(yHTCyL)

# Hair-trap sampling Cantabria
bearHTCA.ch <- read.capthist("./data/captHairTrapsCA.txt", "./data/trapPeloCA.txt", detector='count', noccasions=4)
bearHTCA<-aperm(bearHTCA.ch,c(1,3,2))
yHTCA<-apply(bearHTCA,c(1,2),sum)
(nind<-dim(yHTCA)[1])
yHTCA[,1]<-0
sum(yHTCA)

sex<-read.table("./data/SEX.txt", header=TRUE)
sex<-as.numeric(as.factor(sex$Sex))
head(sex)

# Effort
EF<-read.table("./data/Effort.txt", header=FALSE)[,1]
Eff<-EF/10000

# Region
CA<-read.table("./data/CCAA.txt", header=TRUE)
CA<-as.numeric(as.factor(CA[,2]))

bearPoly<-vect('./GIS/bearPoly.shp')
elev<-rast('./GIS/elevation')
elev<-aggregate(elev,8)
names(elev)<-"Elev"
elev <- crop(elev, bearPoly, mask=TRUE)
elevC<-elev

library(raster)
# X
rx <- raster(ncol=405, nrow=253, xmn=-313100.2, xmx=1306900, ymn=3942375, ymx=4954375)
values(rx)<-matrix(rep(1:405,253),ncol=405,nrow=253, byrow=TRUE)
crs(rx)<-crs(elev)
rx<-rast(rx)
Xc<- crop(rx, bearPoly, mask=TRUE)
# Xc<-raster(Xc)

# Y
ry <- raster(ncol=405, nrow=253, xmn=-313100.2, xmx=1306900, ymn=3942375, ymx=4954375)
values(ry)<-matrix(rep(1:253,405),ncol=405,nrow=253, byrow=FALSE)
crs(ry)<-crs(elev)
ry<-rast(ry)
Yc<- crop(ry, bearPoly, mask=TRUE)
# Yc<-raster(Yc)

elev<-raster(elev)
# Centroids
traplocs<-secr::traps(bearS.ch)

# Castilla y León hair-traps
traplocs2<-traps(bearHTCyL.ch)
Xht1<-data.matrix(traplocs2)
rownames(Xht1)<-1:706
colnames(Xht1)<-c("X","Y")

# Cantabria hair-traps
traplocs3<-traps(bearHTCA.ch)
Xht2<-data.matrix(traplocs3)
rownames(Xht2)<-1:45
colnames(Xht2)<-c("X","Y")

traps<-rbind(data.frame(traplocs),data.frame(traplocs2),data.frame(traplocs3))
plot(bearPoly,asp=TRUE)
points(traps, pch="+")

Xs<-raster(Xc)
Ys<-raster(Yc)
JJ<-stack(elev,Xs,Ys)
names(JJ)<-c('Elev','X','Y')

mymask <- convertRaster(JJ, traps)
str(mymask)
X<-mymask$trapMat
habMat<-mymask$habMat
area<-mymask$area

Xc<-mymask$X
Xc<-rast(Xc)
plane_coord<-as.data.frame(Xc, xy=TRUE)

nsite<-80*44
# The temporary GAM we will take apart.
tmp_jags <- mgcv::jagam(
  response ~ s(x,k=10)+s(y,k=5),
  data = data.frame(
    response = rep(1, nsite),
    x = plane_coord[,1],
    y = plane_coord[,2]),
  family = "poisson",
  file = "tmp.jags"
)

str(tmp_jags$jags.data$S1)
str(tmp_jags$jags.data$X)

elev<-rast(elev)

XX<-array(tmp_jags$jags.data$X,c(44,80,15))

# Create raster layers D1 to D9 using a loop
D_list <- list()
for (i in 2:14) {
  D <- rast(t(XX[,,i]))                     # Transpose and convert to SpatRaster
  crs(D) <- "epsg:25830"                    # Set coordinate reference system
  ext(D) <- ext(elevC)                      # Match extent to elevation raster
  D <- resample(D, elevC)                   # Resample to match resolution
  D <- crop(D, elev, mask = TRUE)           # Crop and mask to study area
  D_list[[i - 1]] <- raster(D)              # Convert to RasterLayer and store
}

# Stack all raster layers into a single object
JJ <- stack(D_list)
names(JJ) <- paste0("D", 1:13)              # Name the layers D1 to D13
JJ <- scale(JJ)                             # Standardize the values
par(mfrow = c(5, 5),mar=c(3.5,3.5,3.5,3.5))  
plot(JJ, ask=TRUE)
dev.off()

mymask <- convertRaster(JJ, as.data.frame(traps))
str(mymask)
dev.off()

## define the model
code <- nimbleCode({

  psi.sex ~ dunif(0,1)
  psi ~ dunif(0, 1)
  beta[1] ~ dnorm(0,0.75)   
  K1[1:9,1:9] <- S1[1:9,1:9] * lambdaS[1]  + S1[1:9,10:18] * lambdaS[2]   
  beta[ 2:10] ~ dmnorm(zeroS[ 2:10],K1[1:9,1:9])
  K2[1:4,1:4] <- S2[1:4,1:4] * lambdaS[3]  + S2[1:4, 5: 8] * lambdaS[4]
  beta[11:14] ~ dmnorm(zeroS[11:14],K2[1:4,1:4]) 
  ## smoothing parameter priors...
  for(i in 1:4){
    lambdaS[i] ~ dgamma(.05,.005)
  } 
  
  alpha1.pS ~ dunif(-10, 10)
  alpha1.pH ~ dunif(-10, 10)

  for (r in 1:4) {
    alpha2.pS[r] ~ dunif(-10, 10)
    alpha2.pH[r] ~ dunif(-10, 10)
  }  
  lp0S.sex ~ dnorm(0, 0.01)
  lp0H.sex ~ dnorm(0, 0.01)
  
  for(i in 1:2){
    p0HT1[i] ~ dunif(0,5)
  }	
  p0HT2 ~ dunif(0,5) 

  # Mixture proportions for males and females
  for(s in 1:2){
    pi_group[s, 1:2] ~ ddirch(dirich_alpha[s, 1:2])  # s = 1 (females), s = 2 (males)
  }

  for(s in 1:2){
    for(g in 1:2){
      alpha1_mix[s,g] ~ dunif(0,10)
      sigma_mix[s,g] <- sqrt(1/(2*alpha1_mix[s,g]))
      sigmaR_mix[s,g] <- sigma_mix[s,g] * pixelWidth
    }
  }
  
  # Constraint: ensure sigma_mix[ ,1] < sigma_mix[ ,2] for both sexes
  one[1] ~ dconstraint(sigma_mix[1,1] < sigma_mix[1,2])  # females
  one[2] ~ dconstraint(sigma_mix[2,1] < sigma_mix[2,2])  # males

  for(j in 1:nTraps1){
    log(p0S[1,j]) <-  alpha1.pS*log.Eff[j] + alpha2.pS[CA[j]] + lp0S.sex
    log(p0S[2,j]) <-  alpha1.pS*log.Eff[j] + alpha2.pS[CA[j]] - lp0S.sex
    log(p0H[1,j]) <-  alpha1.pH*log.Eff[j] + alpha2.pH[CA[j]] + lp0H.sex
    log(p0H[2,j]) <-  alpha1.pH*log.Eff[j] + alpha2.pH[CA[j]] - lp0H.sex
  }

  # Spatial covariates for density
  for(i in 1:(upperLimit[1]-1)) {
    for(j in 1:(upperLimit[2]-1)) {
      log(lam[i, j]) <- inprod(beta[2:10], XG[i,j,1:9]) + inprod(beta[11:14], XG[i,j,10:13])
    }
  }

  lamZ[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))] <- lam[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  * habMat[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]        # convert 'lam' to 0 for non-habitat
  probs[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  <- lamZ[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  / sum(lamZ[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))] ) # 'probs' must sum to 1
  
  for (i in 1:M){
    SEX[i] ~ dbern(psi.sex)
    SEX2[i]<-SEX[i] + 1
    group[i] ~ dcat(pi_group[SEX2[i], 1:2])
    alpha1_ind[i] <- alpha1_mix[SEX2[i], group[i]]
    sigma_ind[i] <- sqrt(1/(2*alpha1_ind[i]))
    males[i] <- z[i] * SEX[i]	        #  Realized males
    females[i] <- z[i] * (1-SEX[i])     #  Realized females
    z[i] ~ dbern(psi)
    S[i, 1] ~ dunif(1, upperLimit[1]) # priors for the activity centres for each individual
    S[i, 2] ~ dunif(1, upperLimit[2])
    negLogDen[i] <- -log(probs[trunc(S[i,1]), trunc(S[i,2])]) # zeros trick
    zerosL[i] ~ dpois(negLogDen[i])
	
    pS[i,1:nTraps1] <- GetDetectionRate1(s = S[i,1:2], 
                                         X = trapMat1[1:nTraps1,1:2], 
                                         J=nTraps1,
                                         alpha1=alpha1_ind[i],
                                         p0=p0S[SEX2[i],1:nTraps1], 
                                         z=z[i])
    pH[i,1:nTraps1] <- GetDetectionRate1(s = S[i,1:2], 
                                         X = trapMat1[1:nTraps1,1:2], 
                                         J=nTraps1,
                                         alpha1=alpha1_ind[i],
                                         p0=p0H[SEX2[i],1:nTraps1], 
                                         z=z[i])
    pHT1[i,1:nTraps2]<-GetDetectionRate2(s = S[i,1:2], 
                                         X = trapMat2[1:nTraps2,1:2], 
                                         J=nTraps2,
                                         alpha1=alpha1_ind[i],
                                         p0=p0HT1[SEX2[i]],
                                         z=z[i])											
    pHT2[i,1:nTraps3]<-GetDetectionRate2(s = S[i,1:2], 
                                         X = trapMat3[1:nTraps3,1:2], 
                                         J=nTraps3,
                                         alpha1=alpha1_ind[i],
                                         p0=p0HT2, 
                                         z=z[i])	
  }
  for(i in 1:n0){
    # Loop through the centroid locations
    # Scat samples
    yS[i,1:nTraps1] ~ dPoissonVector(pS[i,1:nTraps1])
    # Hair samples
    yH[i,1:nTraps1] ~ dPoissonVector(pH[i,1:nTraps1])
    # Hair-traps samples in Castilla y León hair traps
    yHT1[i,1:nTraps2] ~ dPoissonVector(pHT1[i,1:nTraps2]) 
    # Hair-traps samples in Cantabria hair traps
    yHT2[i,1:nTraps3] ~ dPoissonVector(pHT2[i,1:nTraps3])
  }
  ## Model for augmented data (R. Chandler)
  for(i in (n0+1):M) {
    zeros1[i] ~ dpois(sum(pS[i,1:nTraps1]))
    zeros2[i] ~ dpois(sum(pH[i,1:nTraps1]))
    zeros3[i] ~ dpois(sum(pHT1[i,1:nTraps2]))
    zeros4[i] ~ dpois(sum(pHT2[i,1:nTraps3]))
  }
  
  N <- sum(z[1:M])               # Realized number of individuals
  Nmales <- sum(males[1:M]) 	 # Realized number of males
  Nfemales <- sum(females[1:M])  # Realized number of females
  SR <- Nmales / N			     # Male sex ratio
  D <- 1e+8*N/area               # Density (individuals/100 sq.km)
})


GetDetectionRate1 <- nimbleFunction(
  run = function(s = double(1), p0=double(1), alpha1=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- p0*exp(-alpha1*d2)
     return(ans)
    }
  }
)
GetDetectionRate2 <- nimbleFunction(
  run = function(s = double(1), p0=double(0), alpha1=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- p0*exp(-alpha1*d2)
     return(ans)
    }
  }
)

# Vectorized Poisson
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))



n0<- nind 

# Polygons centroids ('traps')
X1<-data.matrix(traplocs)
rownames(X1)<-1:668
colnames(X1)<-c("X","Y")

# Define data augmentation
M<-500
SEX <- c(sex - 1, rep(NA, (M-nind)))

YaugS<-array(0,c(M,nrow(X1)))
YaugS[1:nind,]<-yS
YaugH<-array(0,c(M,nrow(X1)))
YaugH[1:nind,]<-yH
YaugHTCyL<-array(0,c(M,nrow(Xht1))); YaugHTCyL[1:nind,]<-yHTCyL
YaugHTCA <-array(0,c(M,nrow(Xht2))); YaugHTCA[1:nind,] <-yHTCA


XG<-array(NA,c(80,44,13))
for (i in 1:13) {
  XG[,,i] <- mymask[[paste0("D", i)]] * habMat
}

# Organise data for Nimble
zeros<-c(rep(NA, n0), rep(0, M-n0))
str(data  <-   list(yS = YaugS,
                    yH = YaugH,
                    yHT1 = YaugHTCyL,
                    yHT2 = YaugHTCA,
                    log.Eff=log(Eff),
                    XG=XG,
                    SEX=SEX,
                    one=c(1,1),
                    S1 = tmp_jags$jags.data$S1,
                    S2 = tmp_jags$jags.data$S2,
                    zeroS = tmp_jags$jags.data$zero,
                    zerosL = rep(0, M),
                    zeros1=zeros,
                    zeros2=zeros,
                    zeros3=zeros,
                    zeros4=zeros,					
                    habMat=mymask$habMat, 
                    trapMat1=mymask$trapMat[1:668,],
                    trapMat2=mymask$trapMat[669:1374,],
                    trapMat3=mymask$trapMat[1375:1419,]))

str(constants<-list(M = M,
                    n0 = n0,
                    CA=CA,
                    dirich_alpha = matrix(rep(1, 4), nrow = 2, ncol = 2),
                    nTraps1 = nrow(mymask$trapMat[1:668,]),
                    nTraps2 = nrow(mymask$trapMat[669:1374,]),
                    nTraps3 = nrow(mymask$trapMat[1375:1419,]),
                    pixelWidth=mymask$pixelWidth, 
                    area=area,
                    upperLimit=mymask$upperLimit))


nTraps<-nrow(X)
ySi<-yHi<-yHTi<-yHT2i<-array(0,c(nind,nTraps))
ySi[1:nind,1:668]<-yS
yHi[1:nind,1:668]<-yH
yHTi[1:nind,669:1374]<-yHTCyL
yHT2i[1:nind,1375:1419]<- yHTCA

yT<-ySi+yHi+yHTi+yHT2i

Sst <- matrix(NA, nrow(yT), 2)
for(i in 1:nrow(yT)) {
  captTraps <- which(yT[i, ] > 0) # Which traps caught the animal
  captLocs <- mymask$trapMat[captTraps, , drop=FALSE] # Locations of the traps
  Sst[i, ] <- colMeans(captLocs)
  # Check it's in good habitat (might not be if trap is on edge):
  stopifnot(mymask$habMat[Sst[i, , drop=FALSE]] == 1)
}
Sst <- randomPoints(M, mymask, Sst)
                                         
# Initial values
set.seed(1960)
zst<-c(rep(1,nind),rbinom((M-nind),1,0.7))
zeros.init<-c(rep(0, n0), rep(NA, M-n0))
alpha1_mix <- matrix(NA, nrow = 2, ncol = 2)
alpha1_mix[1, ] <- sort(runif(2, 0.05, 2), decreasing = TRUE)  # females: alpha1[1] > alpha1[2] → sigma[1] < sigma[2]
alpha1_mix[2, ] <- sort(runif(2, 0.05, 2), decreasing = TRUE)  # males: alpha1[1] > alpha1[2] → sigma[1] < sigma[2]
str(inits  <-  list(z = zst, 
                    S = Sst,
                    beta= tmp_jags$jags.ini$b,
                    lambdaS = tmp_jags$jags.ini$lambda,
                    alpha1.pS=runif(1, 0, 2),
                    alpha1.pH=runif(1, 0, 2),
                    alpha2.pS=runif(4, 0, 2),
                    alpha2.pH=runif(4, 0, 2),
                    lp0S.sex = runif(1, -2, 2),
                    lp0H.sex = runif(1, -2, 2),
                    p0HT1=runif(2,0,0.2),
                    p0HT2=runif(1,0,0.2),
                    pi_group = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE),
                    alpha1_mix = alpha1_mix,
                    group = sample(1:2, M, replace = TRUE),					
                    SEX=c(rep(NA, nind), rbinom((M-nind), 1, 0.5)),
                    zeros1=zeros.init,
                    zeros2=zeros.init,
                    zeros3=zeros.init,
                    zeros4=zeros.init,
                    psi.sex=runif(1,0.4,0.7),
                    psi=runif(1,0.6,0.9)))

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check=FALSE, calculate=FALSE)
# Function to initialize complex nodes
InitNod<-function(simNodes){
  simNodeScalar <- Rmodel$expandNodeNames(simNodes)
  allNodes <- Rmodel$getNodeNames()
  nodesSorted <- allNodes[allNodes %in% simNodeScalar]
  set.seed(1960) # to fix simulations
  for(n in nodesSorted) {
    Rmodel$simulate(n)
    depNodes <- Rmodel$getDependencies(n)
    Rmodel$calculate(depNodes)
  }
}

InitNod(simNodes = 'SEX2')
any(is.na(Rmodel$SEX2))

Rmodel$initializeInfo()
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
params<-c('N','Nmales','Nfemales','D','SR', 'psi.sex', 
          'alpha1.pS','alpha1.pH',
          'alpha2.pS','alpha2.pH',
          'lp0S.sex','lp0H.sex',
          'sigmaR_mix', 'psi',
          'beta')
params2<-c('S','z','group','SEX')
		  
conf<-configureMCMC(Rmodel, monitors=params, monitors2=params2, useConjugacy=FALSE, enableWAIC = TRUE)

# Rebuild and compile with new sampler
conf$removeSamplers("S")
ACnodes <- paste0("S[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  conf$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}
conf$removeSampler("beta[2:10]")
conf$addSampler(target = paste0("beta[2:10]"), type = "RW_block")
conf$removeSampler("beta[11:14]")
conf$addSampler(target = paste0("beta[11:14]"), type = "RW_block")
conf$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) conf$addSampler(target = node, type = 'slice')
conf$removeSamplers("alpha1_mix")
for(node in Rmodel$expandNodeNames("alpha1_mix")) conf$addSampler(target = node, type = "slice")

MCMC <- buildMCMC(conf)
CompMCMC <- compileNimble(MCMC, project = Rmodel)


## Execute MCMC algorithm and extract samples
## Run model
nb=50000       # Burnin
ni=50000 +nb   # Iters
nc=3          # Chains

start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                  setSeed = TRUE, progressBar = TRUE, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # post-compilation run time
save(outNim, file="Oso_UNIF.RData")


mc<-mcmcOutput(window(outNim$samples,1))
mc2<-mcmcOutput(window(outNim$samples2,1))
summary(mc)
diagPlot(mc)

# MCMC values from mcmc object ‘outNim1$samples’ 
# The object has 37 nodes with 10000 draws for each of 1 chains.
# l95 and u95 are the limits of a 95% Highest Density Credible Interval.
# MCEpc is the Monte Carlo standard error as a percentage of the posterior SD:
        # largest is 9.9%; 21 (57%) are greater than 5.

                     # mean      sd   median      l95      u95 MCEpc
# D                   1.090   0.059    1.085    0.977    1.205 6.896
# N                 371.714  20.243  370.000  333.000  411.000 6.896
# Nfemales          167.054  15.511  166.000  137.000  196.000 6.061
# Nmales            204.661  10.363  204.000  184.000  224.000 5.524
# SR                  0.551   0.024    0.552    0.505    0.597 4.267
# alpha1.pH           1.001   0.105    1.000    0.807    1.220 3.872
# alpha1.pS           0.888   0.066    0.887    0.764    1.023 4.771
# alpha2.pH[1]       -4.208   0.339   -4.211   -4.850   -3.548 3.085
# alpha2.pH[2]       -3.109   0.444   -3.085   -3.999   -2.270 2.589
# alpha2.pH[3]       -2.353   0.155   -2.346   -2.649   -2.052 4.584
# alpha2.pH[4]       -3.979   1.312   -3.784   -6.690   -1.797 2.514
# alpha2.pS[1]       -1.232   0.135   -1.231   -1.492   -0.970 4.966
# alpha2.pS[2]       -0.025   0.161   -0.026   -0.341    0.291 4.702
# alpha2.pS[3]       -0.871   0.091   -0.870   -1.052   -0.699 4.779
# alpha2.pS[4]       -1.687   0.627   -1.660   -3.024   -0.515 3.823
# beta[1]            -0.010   1.162    0.002   -2.255    2.219 0.940
# beta[2]            -1.034   0.779   -1.065   -2.634    0.338 9.436
# beta[3]             1.358   2.282    0.954   -1.987    5.926 9.915
# beta[4]            -0.897   0.832   -0.673   -2.770    0.215 9.706
# beta[5]             1.475   1.226    1.620   -0.950    3.554 9.864
# beta[6]            -0.389   0.853   -0.524   -1.899    1.412 9.578
# beta[7]            -0.482   1.055   -0.630   -2.164    1.504 9.838
# beta[8]             0.206   0.862   -0.003   -1.026    2.098 9.684
# beta[9]             0.059   3.966   -0.538   -5.902    7.941 9.926
# beta[10]            0.092   0.754   -0.015   -1.235    2.084 9.166
# beta[11]            1.461   1.423    1.103   -0.566    5.095 8.750
# beta[12]           -1.467   1.466   -1.401   -4.543    1.308 7.969
# beta[13]            3.712   1.833    3.674   -0.070    7.339 8.072
# beta[14]           -0.815   1.589   -0.216   -5.142    1.239 8.726
# lp0H.sex           -0.217   0.122   -0.215   -0.455    0.016 3.854
# lp0S.sex            0.365   0.068    0.365    0.232    0.504 4.383
# psi                 0.742   0.045    0.740    0.655    0.827 6.515
# psi.sex             0.551   0.036    0.551    0.484    0.622 3.860
# sigmaR_mix[1, 1] 2312.800 144.761 2303.857 2058.098 2626.372 5.332
# sigmaR_mix[2, 1] 4125.206 212.228 4122.050 3715.201 4526.565 6.225
# sigmaR_mix[1, 2] 5111.677 341.034 5087.571 4460.073 5769.446 4.167
# sigmaR_mix[2, 2] 8530.219 527.623 8481.947 7570.037 9584.095 5.985

