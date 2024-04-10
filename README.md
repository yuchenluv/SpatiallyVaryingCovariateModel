# Spatially Varying Covariate Model

A detailed description of the scripts
In scripts

get_data.jl
Functions to process AMS precipitation and climate covariates

params.jl
Constraints for extracting AMS precipitation

generate_data.qmd
Transfer data to be csv file, and as input to the R stan model.

Nonpooled-Nonstationary-Model
Nonstationary model for individual stations. Stan model and R code to run the model.

Pooled-Stationary-Model
Spatially pooled stationary model. Stan model and R code to run the model.

Spatially-Varying-Covariate-Model
The final model. Stan model and R code to run the model.

util.jl
Utility functions for processing model outputs and plotting.

plot_maps.qmd
Plot model outputs.

validation_plots.qmd
Plots for validation.