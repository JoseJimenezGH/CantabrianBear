# ================================================================
# Variation Partitioning of Expected Intensity from SCR Model
# ================================================================
# This script performs variation partitioning analysis on the expected intensity 
# (λ[i,j]) derived from a Spatial Capture-Recapture (SCR) model. It quantifies 
# the relative contribution of covariate groups to spatial variation in predicted intensity.
#
# Key components:
# - Loading posterior samples and covariate matrices
# - Computing posterior means of expected intensity (λ[i,j])
# - Structuring covariates into thematic groups: GAM, Elevation, Forest, Interaction
# - Applying variation partitioning using the 'vegan' package
# - Visualizing and summarizing the results
#
# Requirements:
# - File: "outNim.RData", "BearData.RData"
# - R packages: coda, vegan, mcmcOutput
#
# Date: October 2025

library(coda)
library(vegan)
library(mcmcOutput)

# Load the data
load("C:/.../outNim.RData")
load("C:/.../BearData.RData")

# Extract covariate matrices
XG <- data$XG
Elev <- data$Elev
Elev2 <- data$Elev2
Forest <- data$Forest
Forest2 <- data$Forest2

# Step 1: Extract expected intensity (lam[i,j]) from the SCR model
samples <- mcmcOutput(outNim$samples)

# Extract parameter names corresponding to lam[i,j]
lam_names <- grep("^lam\\[", colnames(samples), value = TRUE)

# Step 2: Compute posterior means of lam[i,j]
# Compute posterior means for beta and b parameters
beta_means <- colMeans(as.matrix(samples[, grep("^beta\\[", colnames(samples))]))
b_means <- colMeans(as.matrix(samples[, grep("^b\\[", colnames(samples))]))

# Convert to spatial matrix
upperLimit <- c(81, 45)
lam_matrix <- matrix(NA, nrow = upperLimit[1]-1, ncol = upperLimit[2]-1)

# Loop through each cell to compute expected intensity
for (i in 1:(upperLimit[1]-1)) {
  for (j in 1:(upperLimit[2]-1)) {
    xb1 <- sum(beta_means[2:10] * XG[i,j,1:9])       # GAM group A
    xb2 <- sum(beta_means[11:14] * XG[i,j,10:13])    # GAM group B
    elev_term <- b_means[1]*Elev[i,j] + b_means[2]*Elev2[i,j]
    forest_term <- b_means[3]*Forest[i,j] + b_means[4]*Forest2[i,j]
    interaction_term <- b_means[5]*Elev[i,j]*Forest[i,j]
    
    log_lam_ij <- xb1 + xb2 + elev_term + forest_term + interaction_term
    lam_matrix[i,j] <- exp(log_lam_ij)
  }
}

# Step 3: Prepare covariate groups
# Response variable: log of expected intensity
n_cells <- (upperLimit[1]-1)*(upperLimit[2]-1)
log_lam <- log(as.vector(lam_matrix))

# GAM group (combined A and B)
GAM <- matrix(NA, nrow = n_cells, ncol = 13)
for (k in 1:13) {
  GAM[,k] <- as.vector(XG[,,k])
}

# Elevation group
Elev_group <- data.frame(
  Elev = as.vector(Elev),
  Elev2 = as.vector(Elev2)
)

# Forest group
Forest_group <- data.frame(
  Forest = as.vector(Forest),
  Forest2 = as.vector(Forest2)
)

# Interaction group
Interaccion <- data.frame(
  ElevForest = as.vector(Elev * Forest)
)

# Check structure of variables
str(log_lam)
str(GAM)
str(Elev_group)
str(Forest_group)
str(Interaccion)

# Step 4: Apply variation partitioning
varpart_result <- varpart(log_lam, GAM, Elev_group, Forest_group)

# Plot the results
plot(varpart_result,
     bg = c("skyblue", "lightgreen", "orange"),  # background colors for each circle
     Xnames = c("GAM", "Elevation", "Forest"),   # group labels
     cutoff = 0.01,                              # threshold to display values
     digits = 2,                                 # number of decimals
     cex = 1.2,                                  # text size
     main = "")

# Show summary of results
summary(varpart_result)

