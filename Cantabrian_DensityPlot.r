# ================================================================
# Spatial Visualization of Density Estimates from SCR Models in R
# ================================================================
# This script performs spatial visualization of individual activity centers 
# derived from Spatial Capture-Recapture (SCR) models using posterior samples.
#
# Key components:
# - Loading required libraries and custom SCR functions
# - Extracting posterior samples of activity centers and latent states
# - Estimating density over a spatial grid
# - Reprojecting and masking raster layers to match study area and bear distribution
# - Creating a final map with contextual geographic layers and annotations
#
# Requirements:
# - Files: "SCR_functions.R", "data/mymask.RData", shapefiles in "GIS" folder
# - R packages: terra, ggplot2, tidyterra, nimble, scrbook, etc.
#
# Date: October 2025

# Load required libraries for spatial analysis, visualization, and SCR modelling
library(terra); library(ggplot2); library(tidyterra); library(png); library(grid)
library(nimble); library(mcmcOutput); library(scrbook); library(scales); library(ggspatial)
library(makeJAGSmask); library(scrbook)

# Load custom SCR functions and pre-generated mask object
source("SCR_functions.R")
load("data/mymask.RData")

# Set working directories (adjust if running on a different system)

#========================
# SPATIAL VISUALIZATION
#========================

# Extract posterior samples of activity centers and latent state z
mco <- mcmcOutput(outNim$samples2)
s <- mco$S
Sx <- as.matrix(s[,,1])  # x-coordinates of activity centers
Sy <- as.matrix(s[,,2])  # y-coordinates of activity centers
z <- mco$z               # latent state (1 = detected individual)

# Bundle into a list for density estimation
obj <- list(Sx = Sx, Sy = Sy, z = z)

# Function to estimate SCR density over a grid and visualize it
SCRdensity <- function(obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
                       Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10) {
  # Extract coordinates and latent states
  Sxout <- obj$Sx
  Syout <- obj$Sy
  z <- obj$z
  niter <- nrow(z)  # number of MCMC iterations
  
  # Define grid boundaries
  if (is.null(Xl)) Xl <- min(Sxout, na.rm = TRUE) * 0.999
  if (is.null(Xu)) Xu <- max(Sxout, na.rm = TRUE) * 1.001
  if (is.null(Yl)) Yl <- min(Syout, na.rm = TRUE) * 0.999
  if (is.null(Yu)) Yu <- max(Syout, na.rm = TRUE) * 1.001
  
  # Create grid
  xg <- seq(Xl, Xu, length.out = nx)
  yg <- seq(Yl, Yu, length.out = ny)
  
  Sx_cut <- cut(Sxout[z == 1], breaks = xg)
  Sy_cut <- cut(Syout[z == 1], breaks = yg)
  
  # Compute density per grid cell
  Dn <- table(Sx_cut, Sy_cut) / niter
  area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
  Dn <- (Dn / area) * scaleout
  
  cat("mean: ", mean(Dn, na.rm = TRUE), fill = TRUE)
  
  # Compute cell centers for plotting
  x_centers <- (xg[-1] + xg[-length(xg)]) / 2
  y_centers <- (yg[-1] + yg[-length(yg)]) / 2
  
  # Convert to matrix format
  Dn_matrix <- matrix(as.numeric(Dn), 
                      nrow = length(x_centers), 
                      ncol = length(y_centers))
  
  # Quick visualization
  par(mar = c(3, 3, 3, 6))
  image(x_centers, y_centers, Dn_matrix, col = terrain.colors(ncolors))
  image.scale(Dn_matrix, col = terrain.colors(ncolors))
  box()
  
  # Return grid and density matrix
  grid <- expand.grid(x = x_centers, y = y_centers)
  return(list(grid = grid, Dn = Dn_matrix))
}

# Run density estimation using dimensions from makeJAGSmask object
res <- SCRdensity(obj, 
                  nx = as.numeric(mymask$upperLimit)[1], 
                  ny = as.numeric(mymask$upperLimit)[2],
                  scalein = 100, scaleout = 100, ncolors = 10)

# Convert density table to numeric matrix
Dn_matrix <- matrix(as.numeric(res$Dn), 
                    nrow = nrow(res$Dn), 
                    ncol = ncol(res$Dn),
                    dimnames = dimnames(res$Dn))

# Rotate matrix 270º counterclockwise (transpose + vertical flip)
Dn_matrix_corr <- t(Dn_matrix)
Dn_matrix_corr <- Dn_matrix_corr[nrow(Dn_matrix_corr):1, ]

# Create raster from corrected matrix
r <- rast(Dn_matrix_corr)

# Assign spatial extent and CRS
xg <- res$grid[,1]
yg <- res$grid[,2]
ext(r) <- ext(min(xg), max(xg), min(yg), max(yg))
crs(r) <- "EPSG:25830"

# Plot raster
plot(r)

# Function to reproject raster to original coordinates using JAGSmask
raster_to_original_coords <- function(r, JAGSmask) {
  if (!inherits(JAGSmask, "JAGSmask")) {
    stop("'", deparse(substitute(JAGSmask)), "' is not a valid 'JAGSmask' object.")
  }
  
  pixWidth <- pixelWidth(JAGSmask)
  origin <- attr(JAGSmask, "origin")
  
  nr <- nrow(r)
  nc <- ncol(r)
  
  x_coords <- origin[1] + (0:(nc - 1)) * pixWidth
  y_coords <- origin[2] + (0:(nr - 1)) * pixWidth
  
  ext(r) <- ext(min(x_coords), max(x_coords), min(y_coords), max(y_coords))
  
  return(r)
}

# Apply coordinate correction
r_corr <- raster_to_original_coords(r, mymask)
plot(r_corr)

# Mask raster to bear distribution polygon
polyBear <- vect('GIS/bearPoly.shp')
Bear <- mask(r_corr, mask = polyBear)
plot(Bear)

# Load additional spatial layers for context
trapShape <- vect("GIS/traps.shp")  # trap locations
prov <- vect('GIS/Provincias.shp')       # province boundaries
AreaEstudio <- vect('GIS/StudyArea.shp')  # study area
ex <- ext(AreaEstudio)

# Load bear image for annotation
bear_img <- rasterGrob(readPNG("data/Designer7.png"), interpolate = TRUE)

# Define custom color palette for density visualization
jet_colors <- colorRampPalette(c("white", "goldenrod", "darkgoldenrod", 
                                 "dark olive green", "darkgreen", "black"))

# Final map with ggplot2
dev.new(width = 11.5, height = 6.2)
ggplot() +
  geom_spatraster(data = Bear) +
  scale_fill_gradientn(colours = jet_colors(1000), na.value = "transparent") +
  geom_spatvector(data = prov, col = "black", fill = NA) +
  geom_point(data = trapShape, aes(x = X, y = Y), shape = "+", size = 2) +
  coord_sf(datum = "EPSG:25830", xlim = c(ex[1], ex[2]), ylim = c(ex[3], ex[4]), expand = FALSE) +
  labs(x = "\nEasting", y = "Northing\n", 
       fill = expression(paste(plain("ind/100 km")^2))) +
  annotation_scale(plot_unit = "m") +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_custom(bear_img, xmin = 340000, xmax = 480000, ymin = 4645000, ymax = 4715000) +
  # Province labels for geographic context
  annotate("text", x = 267889.8, y = 4710000, label = "León", fontface = "italic") +
  annotate("text", x = 369092.2, y = 4714211.3, label = "Palencia", fontface = "italic") +
  annotate("text", x = 285520.2, y = 4798014.9, label = "Asturias", fontface = "italic") +
  annotate("text", x = 219267.9, y = 4665680.1, label = "Zamora", fontface = "italic") +
  annotate("text", x = 150384.2, y = 4760111.6, label = "Lugo", fontface = "italic") +
  annotate("text", x = 158171.4, y = 4688980.7, label = "Orense", fontface = "italic") +
  annotate("text", x = 400000, y = 4786403.3, label = "Cantabria", fontface = "italic") +
  annotate("text", x = 412361.8, y = 4713206.9, label = "Burgos", fontface = "italic") +
  theme_minimal(base_size = 12) +
  theme(panel.background = element_rect(fill = "transparent", colour = "#A4A4A4"),
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank(),

        plot.margin = margin(0.75, 0.45, 0.45, 0.45, "cm"))
