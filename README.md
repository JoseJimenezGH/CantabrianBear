# CantabrianBear
Population Estimation of Cantabrian Brown Bears.

*Bear-SCR-GAM.R*
This R script builds and runs a Spatial Capture-Recapture (SCR) model using the NIMBLE framework to estimate brown bear density in the Cantabrian Mountains (northern Spain). It incorporates genotyped data from scat and hair samples collected along transects, as well as hair samples collected at hair traps. Spatial covariates for abundance are derived from raster layers—including elevation, forest cover, and GAM-based smooth surfaces. The model is formulated as a hierarchical Bayesian framework, with sex- and group-specific detection parameters. 

*DensityPlot.R*
This R script visualizes spatial density estimates from a Spatial Capture-Recapture (SCR) model.
It extracts posterior samples of individual activity centers and latent states from MCMC output.
Density is estimated over a spatial grid and converted into a raster layer with original coordinates.

*SCR_bear_functions*
SCR and plot utilities
