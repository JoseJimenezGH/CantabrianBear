# ================================================================
# Integrated SCR Model for Brown Bear Density Estimation in Spain
# ================================================================
# This script implements a Spatial Capture-Recapture (SCR) model using NIMBLE
# to estimate brown bear density across the Cantabrian Mountains, covering the 
# regions of Asturias, Cantabria, Castilla y León, and Galicia. It integrates 
# scat and hair samples collected along transects, as well as hair-trap samples 
# from Castilla y León and Cantabria. Spatial covariates derived from raster layers 
# are incorporated to account for habitat heterogeneity and environmental gradients.
#
# Key components:
# - Loading and preprocessing of detection histories and spatial covariates
# - Construction of habitat mask and spatial grid using `makeJAGSmask`
# - Generation of GAM-based smooth covariates via `mgcv::jagam`
# - Model specification in NIMBLE with hierarchical structure and spatial intensity
# - Custom detection functions and vectorized Poisson likelihood
# - Data augmentation and initialization of latent variables
# - MCMC configuration with custom samplers and posterior diagnostics
#
# Requirements:
# - R packages: nimble, scrbook, terra, raster, secr, jagsUI, mcmcOutput, tidyverse, ggplot2, etc.
# - External files: detection histories, trap locations, spatial layers (elevation, forest), and custom functions
# - Custom functions: SCR_bear_functions.R, convertRaster(), GetDetectionRate1(), GetDetectionRate2(), dPoissonVector
#
# Output:
# - Posterior estimates of density (D), abundance (N), sex ratio (SR), etc.
# - MCMC diagnostics and summaries
#
# Date: October 2025

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

collS<-data.matrix(read.table("./data/colScats.txt", header=TRUE))[,1]
collH<-data.matrix(read.table("./data/colHair.txt", header=TRUE))[,1]

bearPoly<-vect('./GIS/bearPoly.shp')
elev<-rast('./GIS/elevation')
elev<-aggregate(elev,8)
names(elev)<-"Elev"
elev <- crop(elev, bearPoly, mask=TRUE)
elevC<-elev

forest<-rast('./GIS/forest.tif')
names(forest)<-"Forest"
crs(forest)<-crs(elev)
forestC <- crop(forest, bearPoly, mask=TRUE)

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

elev<-raster(elevC)
forest<-raster(forestC)
forest<-resample(forest,elev)
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
JJ<-stack(elev,forest,Xs,Ys)
names(JJ)<-c('Elev','Forest','X','Y')

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


Elev<-JJ$Elev
Elev2<-(JJ$Elev^2)
Forest<-JJ$Forest
Forest2<-(JJ$Forest)^2

D_list<-c(D_list, Elev, Elev2, Forest, Forest2)

# Stack all raster layers into a single object
JJ <- stack(D_list)
JJ
names(JJ)[1:13] <- paste0("D", 1:13)                        # Name the layers D1 to D13
names(JJ)[14:17]<- c('Elev', 'Elev2', 'Forest', 'Forest2')
JJ <- scale(JJ)                                                 # Standardize the values
par(mfrow = c(6, 5),mar=c(3.5,3.5,3.5,3.5))  
plot(JJ)
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
  
  for(i in 1:5){
    b[i] ~ dnorm(0, 0.01)
  }
  alpha1.pS ~ dunif(-10, 10)
  alpha1.pH ~ dunif(-10, 10)

  for (r in 1:4) {
    alpha2.pS[r] ~ dunif(-10, 10)
    alpha2.pH[r] ~ dunif(-10, 10)
  }
  alpha3.p ~ dunif(-10, 10)
  
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
    log(p0S[1,j]) <-  alpha1.pS*log.Eff[j] + alpha2.pS[CA[j]] + alpha3.p*collS[j] + lp0S.sex
    log(p0S[2,j]) <-  alpha1.pS*log.Eff[j] + alpha2.pS[CA[j]] + alpha3.p*collS[j] - lp0S.sex
    log(p0H[1,j]) <-  alpha1.pH*log.Eff[j] + alpha2.pH[CA[j]] + alpha3.p*collH[j] + lp0H.sex
    log(p0H[2,j]) <-  alpha1.pH*log.Eff[j] + alpha2.pH[CA[j]] + alpha3.p*collH[j] - lp0H.sex
  }
   
  
  # Spatial covariates for intensity
  for(i in 1:(upperLimit[1]-1)) {
    for(j in 1:(upperLimit[2]-1)) {
      log(lam[i, j]) <- inprod(beta[2:10], XG[i,j,1:9]) + inprod(beta[11:14], XG[i,j,10:13]) + b[1]*Elev[i,j] + b[2]*Elev2[i,j] + b[3]*Forest[i,j] + b[4]*Forest2[i,j] + b[5]*Elev[i,j]*Forest[i,j]
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
    S[i, 1] ~ dunif(1, upperLimit[1]) # uniform priors for the activity centres for each individual
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
					Elev=mymask$Elev,
					Elev2=mymask$Elev2,
					Forest=mymask$Forest,
					Forest2=mymask$Forest2,
                    collS=collS,
                    collH=collH,
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
					b=runif(5,-2,2),
                    lambdaS = tmp_jags$jags.ini$lambda,
                    alpha1.pS=runif(1, 0, 2),
                    alpha1.pH=runif(1, 0, 2),
                    alpha2.pS=runif(4, 0, 2),
                    alpha2.pH=runif(4, 0, 2),
					alpha3.p =runif(1,-2, 2),
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
		  'alpha3.p',
          'lp0S.sex','lp0H.sex',
          'sigmaR_mix', 'psi',
          'beta','b')
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

nb=50000
ni=150000 +nb
nc=3
start.time2<-Sys.time()
outNim3 <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains = nc, inits=inits,
                   setSeed = 3, progressBar = TRUE, samplesAsCodaMCMC = TRUE, WAIC = TRUE)
end.time<-Sys.time()
end.time-start.time2


mc<-mcmcOutput(window(outNim$samples,1))
mc2<-mcmcOutput(window(outNim$samples2,1))
summary(mc)
diagPlot(mc)


# MCMC values from mcmc.list object ‘outNim$samples’ 
# The object has 43 nodes with 150000 draws for each of 3 chains.
# l95 and u95 are the limits of a 95% Highest Density Credible Interval.
# Rhat is the estimated potential scale reduction factor:
        # largest is 1.02; NONE are greater than 1.10.
# MCEpc is the Monte Carlo standard error as a percentage of the posterior SD:
        # largest is 3.7%; NONE are greater than 5%.

                     # mean      sd   median      l95      u95  Rhat MCEpc
# D                   1.076   0.062    1.071    0.959    1.197 1.006 1.920
# N                 366.768  21.299  365.000  327.000  408.000 1.006 1.920
# Nfemales          163.828  15.052  163.000  136.000  193.000 0.991 1.479
# Nmales            202.940  11.286  202.000  181.000  224.000 1.012 1.766
# SR                  0.554   0.023    0.554    0.509    0.597 1.000 0.734
# alpha1.pH           0.984   0.109    0.984    0.772    1.199 1.000 0.743
# alpha1.pS           0.861   0.070    0.860    0.727    0.999 1.001 0.974
# alpha2.pH[1]       -4.207   0.343   -4.194   -4.880   -3.536 0.999 0.523
# alpha2.pH[2]       -3.109   0.450   -3.081   -4.005   -2.252 1.001 0.481
# alpha2.pH[3]       -2.356   0.155   -2.353   -2.656   -2.049 1.000 0.753
# alpha2.pH[4]       -3.876   1.301   -3.683   -6.552   -1.625 1.000 0.484
# alpha2.pS[1]       -1.227   0.137   -1.227   -1.496   -0.958 1.000 1.002
# alpha2.pS[2]       -0.014   0.159   -0.014   -0.326    0.300 1.000 0.939
# alpha2.pS[3]       -0.867   0.091   -0.866   -1.047   -0.690 1.000 0.912
# alpha2.pS[4]       -1.597   0.626   -1.566   -2.855   -0.402 1.000 0.797
# alpha3.p            1.039   0.365    1.043    0.303    1.739 1.000 0.609
# b[1]                2.782   1.534    2.715   -0.077    5.816 1.016 3.730
# b[2]               -1.678   1.181   -1.637   -4.008    0.541 1.018 3.727
# b[3]                1.218   0.898    1.148   -0.493    2.948 0.999 3.611
# b[4]               -0.664   0.673   -0.614   -1.987    0.598 1.002 3.598
# b[5]               -0.387   0.249   -0.371   -0.864    0.097 1.000 3.275
# beta[1]            -0.002   1.156   -0.001   -2.278    2.252 1.000 0.144
# beta[2]            -0.878   1.076   -0.971   -2.829    1.488 0.994 2.368
# beta[3]             1.650   2.126    1.448   -2.264    6.096 0.996 2.205
# beta[4]            -0.906   0.838   -0.797   -2.592    0.555 1.000 2.126
# beta[5]             1.185   1.152    1.263   -1.145    3.367 0.998 2.148
# beta[6]            -0.057   0.744   -0.113   -1.399    1.514 0.998 1.875
# beta[7]            -0.400   0.987   -0.471   -2.238    1.627 0.999 2.079
# beta[8]             0.405   0.833    0.258   -0.944    2.103 1.001 2.212
# beta[9]             0.380   3.708    0.055   -6.621    7.947 0.998 2.209
# beta[10]           -0.360   2.039   -0.038   -5.465    3.605 0.988 2.692
# beta[11]            0.594   0.651    0.530   -0.591    1.976 1.001 1.709
# beta[12]           -0.550   1.034   -0.481   -2.753    1.501 1.003 1.526
# beta[13]            1.663   1.323    1.533   -0.782    4.509 1.002 1.768
# beta[14]            0.071   0.642    0.052   -1.265    1.498 0.999 1.719
# lp0H.sex           -0.221   0.122   -0.218   -0.459    0.021 1.000 0.632
# lp0S.sex            0.361   0.067    0.362    0.228    0.491 1.000 0.699
# psi                 0.732   0.047    0.730    0.644    0.826 1.002 1.778
# psi.sex             0.554   0.034    0.554    0.487    0.621 1.001 0.656
# sigmaR_mix[1, 1] 2324.420 141.355 2320.386 2049.578 2603.165 1.000 0.894
# sigmaR_mix[2, 1] 4123.643 213.159 4123.569 3706.251 4543.057 1.000 1.253
# sigmaR_mix[1, 2] 5147.179 346.828 5129.434 4487.204 5838.045 1.000 0.739
# sigmaR_mix[2, 2] 8528.733 521.119 8478.803 7567.296 9589.491 1.000 1.150




