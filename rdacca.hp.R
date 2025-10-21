# ================================================================
# Hierarchical Partitioning of Predictors Using rdacca.hp
# ================================================================
# This script performs hierarchical partitioning to assess the relative 
# importance of predictor blocks in explaining variation in expected intensity 
# (log(λ)) derived from a Spatial Capture-Recapture (SCR) model.
#
# Key components:
# - Loading covariate data and posterior estimates
# - Structuring predictors into thematic blocks: GAM, Elevation, Forest
# - Applying hierarchical partitioning with the 'rdacca.hp' package
# - Aggregating and visualizing the relative contributions of each block
#
# Requirements:
# - Files: "BearData.RData", "varpartVegan.R"
# - R packages: rdacca.hp, GFA, dplyr, ggplot2
#
# Author: Jose Jimenez Garcia-Herrera
# Date: October 2025

# Load required libraries
library(rdacca.hp)
library(GFA)
library(dplyr)
library(ggplot2)

# Purpose: Perform hierarchical partitioning to assess the relative importance of each predictor.
# Load the data
load("C:/.../BearData.RData")
source('C:/.../varpartVegan.R')

# 1. Data preparation
# Convert elevation and forest matrices to vectors
elev_vec <- as.vector(Elev)
forest_vec <- as.vector(Forest)

# Convert 3D array XG into a 2D matrix: rows = pixels, columns = GAM variables
XG_mat <- matrix(NA, nrow = dim(XG)[1] * dim(XG)[2], ncol = dim(XG)[3])
for (i in 1:dim(XG)[3]) {
  XG_mat[, i] <- as.vector(XG[,,i])
}
colnames(XG_mat) <- paste0("XG", 1:ncol(XG_mat))

# Build full data frame
datos <- data.frame(
  log_lam = log_lam,
  elev = elev_vec,
  forest = forest_vec
)
datos <- cbind(datos, XG_mat)

# 2. Define predictor blocks
bloques <- list(
  GAM = colnames(XG_mat),
  Elevation = "elev",
  Forest = "forest"
)

# 3. Create predictor matrix
predictors <- datos[, c(bloques$GAM, bloques$Elevation, bloques$Forest)]

# 4. Run rdacca.hp
resultado <- rdacca.hp(
  datos$log_lam,   # response variable
  predictors,      # predictor matrix
  type = "adjR2",  # use adjusted R²
  n.perm = 999     # number of permutations
)

# Extract hierarchical partitioning results
tabla <- as.data.frame(resultado$Hier.part)
tabla$Variable <- rownames(resultado$Hier.part)

# 5. Aggregate importance by block
tabla_bloques <- tabla %>%
  mutate(Block = case_when(
    Variable %in% colnames(XG_mat) ~ "GAM",
    Variable == "elev" ~ "Elevation",
    Variable == "forest" ~ "Forest"
  )) %>%
  group_by(Block) %>%
  summarise(Importance = sum(Individual))  # Use the 'Individual' column for aggregation
# 6. Display results
print(tabla_bloques)
# A tibble: 3 × 2
#   Block     Importance
#   <chr>          <dbl>
# 1 Elevation     0.194 
# 2 Forest        0.0936
# 3 GAM           0.653 

# 7. Plot results
ggplot(tabla_bloques, aes(x = Block, y = Importance, fill = Block)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(Importance, 3)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("skyblue", "lightgreen", "tan")) +
  labs(title = "Relative Importance by Block",
       x = "",
       y = "Adjusted R² Contribution") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
