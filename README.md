# CantabrianBear
Cantabrian bear population size estimate

*Bear-SCR-GAM.R*
This R script builds and runs a Spatial Capture-Recapture (SCR) model using the NIMBLE framework to estimate brown bear density in the Cantabrian Mountains (N Spain). It processes detection data from scat, hair, and hair-trap samples, constructs spatial covariates from raster layers (elevation, forest, GAM smooths), and defines a hierarchical Bayesian model with sex- and group-specific detection parameters. 

*DensityPlot.R*
This R script visualizes spatial density estimates from a Spatial Capture-Recapture (SCR) model.
It extracts posterior samples of individual activity centers and latent states from MCMC output.
Density is estimated over a spatial grid and converted into a raster layer with original coordinates.

*VarPartVegan.R*
This script quantifies the contribution of different covariate groups (GAM, Elevation, Forest) to spatial variation in expected intensity (λ) from an SCR model. It computes posterior means of model parameters and structures covariates into thematic blocks.

*rdacca.hp.R*
This script applies hierarchical partitioning to quantify the relative importance of predictor blocks (GAM, Elevation, Forest) in explaining expected intensity (log(λ)) from an SCR model. Covariates are grouped and analyzed using the rdacca.hp package with permutation-based adjusted R². Results are aggregated by block and visualized with a bar plot showing their contributions.

*Counterfactual.R*
This script estimates the relative contribution of GAM, Elevation, and Forest covariates to expected intensity (λ) in an SCR model. It simulates counterfactual scenarios by excluding each covariate group from the model and quantifies the reduction in predicted intensity. Results are summarized with uncertainty intervals and visualized using boxplots.

*SCR_bear_functions*
SCR and plot utilities
