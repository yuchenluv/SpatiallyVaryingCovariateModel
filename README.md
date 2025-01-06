# Spatially Varying Covariate Model

## 
To reproduce the results



## 
A detailed description of the scripts
In scripts

get_data.jl
Functions to process AMS precipitation (from GHCN dataset) and climate covariates (CO2 as an example)

params.jl
Constraints for extracting AMS precipitation

get_Altas14.jl
Web scrape the Atlas 14 IDF estimates for comparison.

cal_util.jl
Utility functions to process estimated distribution parameters to return level estimates.

plot_util.jl
Utility functions for plotting.

generate_data.qmd
Transfer data to be csv file, and as input to the R stan model.

Nonpooled-Nonstationary-Model
Nonstationary model for individual stations. Stan model and R code to run the model.

Pooled-Stationary-Model
Spatially pooled stationary model. Stan model and R code to run the model.

Spatially-Varying-Covariate-Model
The final model. Stan model and R code to run the model.

validation_metrics.qmd
Different validation metrics.

validation_plots.qmd
Validation plots including comparison with Atlas 14, cross-validation among different subsets, comparing the three frameworks.