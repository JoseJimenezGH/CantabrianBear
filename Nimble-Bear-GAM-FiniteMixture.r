#==============================================================================#
#                                                                              #
#                            CANTABRIC BROWN BEAR                              #
#                            ~~~ Ursus arctos ~~~                              #
#                               José Jiménez                                   #
#                                 IREC-CSIC                                    #
#                             11:23 14/12/2024                                 #
#                                                                              #
#==============================================================================#

# Load packages
library(nimble)
library(coda)
library(raster)
library(secr)
library(jagsUI)
library(makeJAGSmask)
library(terra)
library(scrbook)

setwd('C:/...')


# SETTING DATA FOR SCR-INTEGRATION
#==================================
# We construct the capture matrices by individual-trap. To do this, for each method we use dummy 
# captures for all uncaptured individuals, which we then remove in a next step. What we have done 
# is to assign the non-captured individuals to a trap with no real captures (e.g. trap '1' in the 
# case of scat) and once we have created the method-matched capture matrices with the same individual 
# order, we have removed the dummy captures (e.g. in the case of scat, by deleting all captures in 
# trap '1'). Similarly for hair traps

# Scat sampling
bearS.ch <- read.capthist("./data/captScat.txt", "./data/traps.txt", detector='count', noccasions=1)
bearS<-aperm(bearS.ch,c(1,3,2))
yS<-apply(bearS,c(1,2),sum)
(nind<-dim(yS)[1])
yS[,1]<-0
sum(yS)
rownames(yS)<-1:nind

# Hair sampling
bearH.ch <- read.capthist("./data/captHair.txt", "./data/traps.txt", detector='count', noccasions=1)
bearH<-aperm(bearH.ch,c(1,3,2))
yH<-apply(bearH,c(1,2),sum)
(nind<-dim(yH)[1])
yH[,1]<-0
sum(yH)

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
sex[,2]<-as.numeric(as.factor(sex$Sex))
head(sex)

## F: 1
## M: 2

# Define data augmentation
M<-500

# Raster data
load("./GIS/allRasters.RData")

# Polygons centroids ('traps')
traplocs<-traps(bearS.ch)
X<-data.matrix(traplocs)
rownames(X)<-1:668
colnames(X)<-c("X","Y")

# We scaled all 'traps' coordinates
X<-X/10000
X[,1]<-X[,1]- minX0
X[,2]<-X[,2]- minY0

# Castilla y León hair-traps
traplocs2<-traps(bearHTCyL.ch)
Xht1<-data.matrix(traplocs2)
rownames(Xht1)<-1:706
colnames(Xht1)<-c("X","Y")
X2<-Xht1/10000
X2[,1]<-X2[,1]- minX0
X2[,2]<-X2[,2]- minY0

# Cantabria hair-traps
traplocs3<-traps(bearHTCA.ch)
Xht2<-data.matrix(traplocs3)
rownames(Xht2)<-1:45
colnames(Xht2)<-c("X","Y")
X3<-Xht2/10000
X3[,1]<-X3[,1]- minX0
X3[,2]<-X3[,2]- minY0

# Data augmentation
YaugS<-array(0,c(M,nrow(X))); YaugS[1:nind,]<-yS
YaugH<-array(0,c(M,nrow(X))); YaugH[1:nind,]<-yH
YaugHTCyL<-array(0,c(M,nrow(X2))); YaugHTCyL[1:nind,]<-yHTCyL
YaugHTCA<-array(0,c(M,nrow(X3))); YaugHTCA[1:nind,]<-yHTCA

table(YaugHTCyL)
table(YaugHTCA)

# We make a data.frame
elev <- as.data.frame(elev.r, xy=TRUE)
gcoords <- elev[,c("x", "y")]
nPix <- nrow(gcoords)
pixelArea <- prod(res(elev.r))

# We stacked and scaled rasters
JJ<-stack(elev.r,Xc,Yc,FOREST)
names(JJ)<-c('Elev','Xc','Yc','FOREST')

par(oma=c(2,2,2,2))
plot(JJ)
JJ<-scale(JJ)

# Region
CA<-read.table("./data/CCAA.txt", header=TRUE)
CA<-as.numeric(as.factor(CA[,2]))

# 1: Asturias
# 2: Cantabria
# 3: Castilla y León
# 4: Galicia

# Effort
EF<-read.table("./data/Effort.txt", header=FALSE)[,1]
Eff<-(EF-mean(EF))/sd(EF)

# Binary variable to indicate whether the sample was collected after a bear attack on an apiary
collS<-data.matrix(read.table("./data/colScats.txt", header=TRUE))[,1]
collH<-data.matrix(read.table("./data/colHair.txt", header=TRUE))[,1]

# We use the habitat mask to limit the analysis to the area where there are bears
mymask <- convertRaster(JJ, as.data.frame(rbind(X,X2,X3)))
str(mymask)
habMat<-mymask$habMat
area<-mymask$area

# load nimble
library(nimble)
## define the model
code <- nimbleCode({

  psi.sex ~ dunif(0,1)
  psi ~ dunif(0, 1)
  
  betaG[1] ~ dnorm(0,0.75)
  ## prior for s(x,y)...
  K1[1:5,1:5] <- S1[1:5, 1:5]* lambdaS[1] + S1[1:5, 6:10]* lambdaS[2]
  betaG[2:6] ~ dmnorm(zeroS[ 2:6],K1[1:5,1:5])
  for(i in 1:2){
    beta[i] ~ dnorm(0,0.01)
  }
  ## smoothing parameter priors...
  for(i in 1:2){
    lambdaS[i] ~ dgamma(.05,.005)
    rho[i] <- log(lambdaS[i])
  } 
  
  alpha1.p ~ dunif(-10, 10)        # Eff (scats)
  alpha2.p ~ dunif(-10, 10)        # Eff^2 (scats)
  lp0S.sex ~ dnorm(0, 0.01)        # RE (sex) scats
  lp0H.sex ~ dnorm(0, 0.01)        # RE (sex) hair
  for(s in 1:2){
    p0HT1[s] ~ dunif(0,1)
  }
  p0HT2 ~ dunif(0,1) 
  
  for (r in 1:nCA) {
    alpha3.pS[r] ~ dunif(-10, 10) # Author (scats)
    alpha3.pH[r] ~ dunif(-10, 10) # Author (hair)
  }
  alpha4.p ~ dunif(-10, 10)
  
  for(s in 1:3){
    sigma[s] ~ dunif(0, 5)
    alpha1[s] <- 1/(2*sigma[s]^2)
    sigmaR[s] <- sigma[s]*pixelWidth
  }
  
  one[1] ~ dconstraint(sigma[1] <= sigma[2])
  one[2] ~ dconstraint(sigma[2] <= sigma[3])
  
  for(j in 1:nTraps1){
    log(p0S[1,j]) <-  alpha1.p*Eff[j] + alpha2.p*(Eff[j])^2  + alpha3.pS[CA[j]] + alpha4.p*collS[j] + lp0S.sex  # Female (scat)
    log(p0S[2,j]) <-  alpha1.p*Eff[j] + alpha2.p*(Eff[j])^2  + alpha3.pS[CA[j]] + alpha4.p*collS[j] - lp0S.sex  # Male (scat)
    log(p0H[1,j]) <-  alpha1.p*Eff[j] + alpha2.p*(Eff[j])^2  + alpha3.pH[CA[j]] + alpha4.p*collH[j] + lp0H.sex  # Female (hair)
    log(p0H[2,j]) <-  alpha1.p*Eff[j] + alpha2.p*(Eff[j])^2  + alpha3.pH[CA[j]] + alpha4.p*collH[j] - lp0H.sex  # Male (hair)
  }

  # Spatial covariates for density
  for(i in 1:(upperLimit[1]-1)) {
    for(j in 1:(upperLimit[2]-1)) {
      log(lam[i, j]) <- inprod(beta[1:2], XG[i,j, 1:2]) + inprod(betaG[2:6], XG[i, j, 3:7])
    }
  }
  lam0[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))] <- lam[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  * habMat[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]       # convert 'lam' to 0 for non-habitat
  probs[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  <- lam0[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))]  / sum(lam0[(1:(upperLimit[1]-1)), (1:(upperLimit[2]-1))] ) # 'probs' must sum to 1

  for (i in 1:M){
    g[i] ~ dcat(pi[1:3]) 
    SEX[i]~dbern(psi.sex)               #  Potential males
    SEX2[i]<-SEX[i] + 1
    males[i] <- z[i] * SEX[i]	        #  Realized males
    sexfemale[i] <- 1-SEX[i]	        #  Potential females
    females[i] <- z[i] * sexfemale[i]   #  Realized females
    z[i] ~ dbern(psi)
    S[i, 1] ~ dunif(1, upperLimit[1]) # uniform priors for the activity centres for each individual
    S[i, 2] ~ dunif(1, upperLimit[2])
    negLogDen[i] <- -log(probs[trunc(S[i,1]), trunc(S[i,2])]) # zeros trick
    zeros[i] ~ dpois(negLogDen[i])
    
	# Detection. Scats in transects
	pS[i,1:nTraps1] <- GetDetectionRate(s = S[i,1:2], 
                                        X = trapMat1[1:nTraps1,1:2], 
                                        J=nTraps1,
                                        alpha1=alpha1[g[i]],
                                        p0=p0S[SEX2[i],1:nTraps1], 
                                        z=z[i])
    # Detection. Hair in transects    
    pH[i,1:nTraps1] <- GetDetectionRate(s = S[i,1:2], 
                                        X = trapMat1[1:nTraps1,1:2], 
                                        J=nTraps1,
                                        alpha1=alpha1[g[i]],
                                        p0=p0H[SEX2[i],1:nTraps1], 
                                        z=z[i])
    # Detection. Hair traps - Castilla y León
	DsqHT1[i,1:nTraps2] <- (S[i,1]-trapMat2[1:nTraps2,1])^2 + (S[i,2]-trapMat2[1:nTraps2,2])^2
    pHT1[i,1:nTraps2] <- p0HT1[SEX2[i]] * exp(-alpha1[g[i]]*DsqHT1[i,1:nTraps2]) * z[i]
	# Detection. Hair traps - Cantabria	
    DsqHT2[i,1:nTraps3] <- (S[i,1]-trapMat3[1:nTraps3,1])^2 + (S[i,2]-trapMat3[1:nTraps3,2])^2	
    pHT2[i,1:nTraps3] <- p0HT2 * exp(-alpha1[g[i]]*DsqHT2[i,1:nTraps3]) * z[i]
    
	# Loop through the centroid locations
    # Scat samples
    yS[i,1:nTraps1] ~ dPoissonVector(pS[i,1:nTraps1])
    # Hair samples
    yH[i,1:nTraps1] ~ dPoissonVector(pH[i,1:nTraps1])
    # Hair-traps samples
    yHT1[i,1:nTraps2] ~ dBinomVector(P=pHT1[i,1:nTraps2], K=4)
    # Hair-traps samples
    yHT2[i,1:nTraps3] ~ dBinomVector(P=pHT2[i,1:nTraps3], K=4)
  }
  N <- sum(z[1:M])               # Realized number of individuals
  Nmales <- sum(males[1:M]) 	 # Realized number of males
  Nfemales <- sum(females[1:M])  # Realized number of females
  SR <- Nmales / N			     # Male sex ratio
  D <- N/area
})


GetDetectionRate <- nimbleFunction(
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


# Vectorized binomial
dBinomVector <- nimbleFunction(
  run = function(x = double(1), P = double(1), K = double(0), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dbinom(x[j], K, P[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rBinomVector  <- nimbleFunction(
  run = function(n = integer(), P = double(1), K = double(0)) {
    J <- length(P)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rbinom(1, K, P[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dBinomVector = list(
    BUGSdist = "dBinomVector(P, K)",
    Rdist = "dBinomVector(P, K)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'P = double(1)', 'K = double(0)'))
))


nTraps1 <- nrow(X)
nTraps2 <- nrow(X2)
nTraps3 <- nrow(X3)

# We prepared the GAM data
Xc<-mymask$Xc
Xc<-raster(Xc)
plane_coord<-as.data.frame(Xc, xy=TRUE)

knots<-6
nsite<-80*45 # raster size

# The temporary GAM we will take apart.
tmp_jags <- mgcv::jagam(
  response ~ y+I(y^2)+s(x,k=knots),
  data = data.frame(
    response = rep(1, nsite),
    x = plane_coord[,1],
    y = plane_coord[,2]),
  family = "poisson",
  file = "tmp.jags"
)

str(tmp_jags$jags.data$S1)
str(tmp_jags$jags.data$X)

load("./GIS/mymask_GAM.RData")

XG<-array(NA,c(80,45,7))
XG[,,1]<-mymask$D1*habMat
XG[,,2]<-mymask$D2*habMat
XG[,,3]<-mymask$D3*habMat
XG[,,4]<-mymask$D4*habMat
XG[,,5]<-mymask$D5*habMat
XG[,,6]<-mymask$D6*habMat
XG[,,7]<-mymask$D7*habMat


SEX <- c(sex[,2] - 1, rep(NA, (M-nind)))

# Organise data for Nimble
str(data  <-   list(yS = YaugS,
                    yH = YaugH,
                    yHT1 = YaugHTCyL,
                    yHT2 = YaugHTCA,
                    zeros = rep(0, M),
                    Eff=Eff, 
                    XG=XG,
                    SEX=SEX,
                    collS=collS,
                    collH=collH,
                    one=c(1,1),
                    S1 = tmp_jags$jags.data$S1,
                    zeroS = tmp_jags$jags.data$zero,
                    habMat=mymask$habMat, 
                    trapMat1=mymask$trapMat[1:668,],
                    trapMat2=mymask$trapMat[669:1374,],
                    trapMat3=mymask$trapMat[1375:1419,]))
str(constants<-list(M = M, 
                    nTraps1 = nTraps1,
                    nTraps2 = nTraps2,
                    nTraps3 = nTraps3,
                    nCA=4,
                    CA=CA,
                    pixelWidth=mymask$pixelWidth, 
                    area=area, 
                    upperLimit=mymask$upperLimit))

# We set up data-compatible inits
nindT<-nrow(yS)
nTraps<-nrow(as.data.frame(rbind(X,X2,X3)))
ySi<-yHi<-yHTi<-yHT2i<-array(0,c(nindT,nTraps))
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
str(inits  <-  list(z = rep(1, M), 
                    S = Sst,
                    beta= tmp_jags$jags.ini$b[2:3], 
                    betaG = tmp_jags$jags.ini$b[1:6],
                    lambdaS = tmp_jags$jags.ini$lambda,
                    sigma=rep(0.3, 3),
                    alpha1.p=runif(1, 0.75, 1.00), 
                    alpha2.p=runif(1,-0.12,-0.05), 
                    alpha3.pS=runif(4,-5.00,-0.50),
                    alpha3.pH=runif(4,-8.80,-3.00),
                    alpha4.p=runif(1, -2, 2), 
                    lp0S.sex =runif(1,-1.00, 0.00),
                    lp0H.sex= runif(1, 0.00, 0.30),
                    p0HT1=runif(2,0,1),
                    p0HT2=runif(1,0,1),
                    pi=rep(1/3,3),
                    g=rep(1,M),
                    SEX=c(rep(NA, nindT), rbinom((M-nindT), 1, 0.5)),
                    psi.sex=runif(1,0.4,0.8),
                    psi=runif(1,0.6,0.9)))

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, calculate=F, check=F)
Rmodel$initializeInfo()
#Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)
params<-c('N','Nmales','Nfemales','D','SR', 'psi.sex', 
          'alpha1.p', 'alpha2.p', 'alpha3.pS', 'alpha3.pH', 'alpha4.p',
          'lp0S.sex', 'lp0H.sex',
          'p0HT1','p0HT2',
          'sigmaR', 'psi', 
          'betaG','beta')
mcmc<-configureMCMC(Rmodel, monitors=params, enableWAIC = TRUE)

# Rebuild and compile with new sampler
mcmc$removeSamplers("S")
ACnodes <- paste0("S[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  mcmc$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

mcmc$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) mcmc$addSampler(target = node, type = 'slice')
 
MCMC <- buildMCMC(mcmc)

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

outNim$WAIC


summary(mcmcOutput(outNim$samples))
diagPlot(mcmcOutput(outNim$samples))


# MCMC values from mcmc.list object ‘outNim$samples’ 
# The object has 34 nodes with 50000 draws for each of 3 chains.
# l95 and u95 are the limits of a 95% Highest Density Credible Interval.
# Rhat is the estimated potential scale reduction factor:
        # largest is 1.03; NONE are greater than 1.10.
# MCEpc is the Monte Carlo standard error as a percentage of the posterior SD:
        # largest is 4.9%; NONE are greater than 5%.

                # mean     sd  median     l95     u95  Rhat MCEpc
# D              1.048  0.055   1.046   0.939   1.152 1.027 2.543
# N            373.869 19.609 373.000 335.000 411.000 1.027 2.543
# Nfemales     138.479  9.256 138.000 121.000 156.000 1.029 1.766
# Nmales       235.389 15.251 235.000 204.000 263.000 1.000 2.407
# SR             0.629  0.019   0.630   0.593   0.667 1.000 1.324
# alpha1.p       0.720  0.063   0.719   0.599   0.846 1.004 2.741
# alpha2.p      -0.074  0.012  -0.074  -0.098  -0.051 1.001 2.601
# alpha3.pH[1]  -3.929  0.354  -3.916  -4.637  -3.254 1.000 1.280
# alpha3.pH[2]  -3.104  0.450  -3.080  -4.006  -2.260 1.000 0.739
# alpha3.pH[3]  -2.314  0.138  -2.313  -2.590  -2.047 1.002 1.455
# alpha3.pH[4]  -3.840  1.301  -3.653  -6.509  -1.602 0.999 0.742
# alpha3.pS[1]  -1.063  0.174  -1.059  -1.411  -0.733 1.000 2.326
# alpha3.pS[2]  -0.029  0.160  -0.030  -0.339   0.287 1.000 1.353
# alpha3.pS[3]  -0.922  0.094  -0.921  -1.112  -0.743 0.999 1.920
# alpha3.pS[4]  -1.593  0.627  -1.565  -2.870  -0.405 1.000 1.053
# alpha4.p       1.102  0.348   1.106   0.406   1.765 1.001 0.878
# beta[1]        0.939  0.136   0.940   0.660   1.192 1.017 4.936
# beta[2]       -0.010  0.001  -0.010  -0.012  -0.007 1.014 4.935
# betaG[1]      -0.001  1.155  -0.006  -2.263   2.273 1.000 0.264
# betaG[2]       1.882  0.890   1.954  -0.187   3.698 0.990 2.413
# betaG[3]      28.306  2.874  28.235  22.787  34.010 1.002 2.273
# betaG[4]       5.204  0.674   5.199   3.845   6.541 1.003 2.281
# betaG[5]      37.964  4.003  37.860  30.220  45.848 0.995 2.266
# betaG[6]       0.219  1.366   0.033  -2.624   3.559 0.996 2.476
# lp0H.sex      -0.361  0.128  -0.358  -0.612  -0.115 0.999 1.537
# lp0S.sex       0.196  0.078   0.199   0.041   0.345 0.999 2.213
# p0HT1[1]       0.002  0.001   0.002   0.001   0.004 0.999 0.887
# p0HT1[2]       0.004  0.001   0.004   0.003   0.006 1.001 0.891
# p0HT2          0.021  0.007   0.020   0.009   0.036 1.001 0.846
# psi            0.747  0.043   0.745   0.662   0.831 1.003 2.358
# psi.sex        0.629  0.031   0.629   0.568   0.689 1.001 1.099
# sigmaR[1]      0.191  0.019   0.190   0.154   0.229 1.000 2.110
# sigmaR[2]      0.308  0.027   0.305   0.261   0.363 0.998 3.047
# sigmaR[3]      0.657  0.021   0.656   0.618   0.699 1.000 2.072
