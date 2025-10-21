# ================================================================
# Counterfactual Reduction Analysis of Expected Intensity (λ)
# ================================================================
# This script evaluates the relative contribution of three covariate groups 
# (GAM, Elevation, Forest) to the expected intensity (λ) in a Spatial 
# Capture-Recapture (SCR) model. It uses posterior samples to simulate 
# counterfactual scenarios where each group is excluded, and quantifies 
# the reduction in predicted intensity.
#
# Key components:
# - Loading posterior samples and covariate matrices
# - Computing expected intensity under full and reduced models
# - Estimating relative reductions and uncertainty intervals
# - Visualizing the distribution of reductions using boxplots
#
# Requirements:
# - File: "outNim.RData", "BearData.RData"
# - R packages: mcmcOutput, ggplot2
#
# Author: Jose Jimenez Garcia-Herrera
# Date: October 2025

library(mcmcOutput)
mc <- mcmcOutput(outNim$samples)

# Load:
load("C:/.../outNim.RData")
load("C:/.../BearData.RData")

# Load the following matrices:
# XG: 3D array [I, J, 13]
# Elev, Elev2, Forest, Forest2: matrices [I, J]
# beta_samples: matrix [n_samples, 14]
# b_samples: matrix [n_samples, 5]

# Posteriors
b_samples<-mc$b
beta_samples<-mc$beta

# Data
XG <- data$XG
Elev <- data$Elev
Elev2 <- data$Elev2
Forest <- data$Forest
Forest2 <- data$Forest2

# Initialize vectors to store results
n_samples <- nrow(beta_samples)
I <- dim(XG)[1]
J <- dim(XG)[2]

lam_total_samples <- numeric(n_samples)
lam_without_GAM_samples <- numeric(n_samples)
lam_without_Elev_samples <- numeric(n_samples)
lam_without_Forest_samples <- numeric(n_samples)

# Loop through posterior samples
for (s in 1:n_samples) {
  beta <- beta_samples[s, ]
  b <- b_samples[s, ]
  
  loglam_smooth1 <- apply(XG[,,1:9], c(1,2), function(x) sum(x * beta[2:10]))
  loglam_smooth2 <- apply(XG[,,10:13], c(1,2), function(x) sum(x * beta[11:14]))
  loglam_elev <- b[1] * Elev + b[2] * Elev2
  loglam_forest <- b[3] * Forest + b[4] * Forest2 + b[5] * Elev * Forest

  lam_total_samples[s] <- mean(exp(loglam_smooth1 + loglam_smooth2 + loglam_elev + loglam_forest))
  lam_without_GAM_samples[s] <- mean(exp(loglam_elev + loglam_forest))
  lam_without_Elev_samples[s] <- mean(exp(loglam_smooth1 + loglam_smooth2 + loglam_forest))
  lam_without_Forest_samples[s] <- mean(exp(loglam_smooth1 + loglam_smooth2 + loglam_elev))
}

# Compute relative reductions
delta_GAM <- mean(lam_total_samples - lam_without_GAM_samples) / mean(lam_total_samples)
delta_Elev <- mean(lam_total_samples - lam_without_Elev_samples) / mean(lam_total_samples)
delta_Forest <- mean(lam_total_samples - lam_without_Forest_samples) / mean(lam_total_samples)

# Display results
data.frame(
  Component = c("GAM", "Elevation", "Forest"),
  Relative_Reduction = c(delta_GAM, delta_Elev, delta_Forest)
)

# Compute reduction vectors for uncertainty analysis
delta_GAM_vec <- (lam_total_samples - lam_without_GAM_samples) / lam_total_samples
delta_Elev_vec <- (lam_total_samples - lam_without_Elev_samples) / lam_total_samples
delta_Forest_vec <- (lam_total_samples - lam_without_Forest_samples) / lam_total_samples

# Function to compute median and 95% CI
summary_ci <- function(x) {
  med <- median(x, na.rm=TRUE)
  ci <- quantile(x, probs = c(0.025, 0.975), na.rm=TRUE)
  return(c(Median = med, CI_2.5 = ci[1], CI_97.5 = ci[2]))
}

# Apply summary function to each component
res_GAM <- summary_ci(delta_GAM_vec)
res_Elev <- summary_ci(delta_Elev_vec)
res_Forest <- summary_ci(delta_Forest_vec)

# Combine results into a data frame
df_summary <- data.frame(
  Component = c("GAM", "Elevation", "Forest"),
  Median = round(c(res_GAM[1], res_Elev[1], res_Forest[1]), 2),
  CI_2.5 = round(c(res_GAM[2], res_Elev[2], res_Forest[2]), 2),
  CI_97.5 = round(c(res_GAM[3], res_Elev[3], res_Forest[3]), 2)
)

# Print summary table
print(df_summary)

# Load ggplot2 for visualization
library(ggplot2)

# Create data frame for plotting
df <- data.frame(
  Component = rep(c("GAM", "Elevation", "Forest"), each = length(delta_GAM_vec)),
  Reduction = c(delta_GAM_vec, delta_Elev_vec, delta_Forest_vec)
)

# Boxplot of counterfactual reductions
ggplot(df, aes(x = Component, y = Reduction, fill = Component)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    title = "Distribution of Counterfactual Reductions",
    y = "Relative Reduction in Expected Intensity",
    x = "Component"
  ) +
  theme_minimal() +
  theme(legend.position = "none")