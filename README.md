# Spatially Varying Covariate Model
This repository includes the codes for the manuscript (...)

## How to use
To reproduce the results:
(1) Preprocess annual maximum rainfall records and CO2 covariate with generate_data.jl, csv files are generated that could be inputs to the statistical model
(2) Choose the framework of interest (Pooled Stationary Model / Nonpooled Nonstationary Model / Spatially Varying Covariate Model), use the corresponding qmd file to estimate the distribution parameters, parameters are saved as csv file as well
(3) One could use plot_results.qmd for most plots
For more details, readers can refer to the following descriptions of the different files

## Detailed descriptions
In scripts

get_data.jl
Functions to process AMS precipitation (from GHCN dataset) and climate covariates (CO2 as an example)

params.jl
Constraints for extracting AMS precipitation

get_Altas14.jl
Web scrape the Atlas 14 IDF estimates for comparison

cal_util.jl
Utility functions to process estimated distribution parameters to return level estimates, also functions for the validation metrics

plot_util.jl
Utility functions for plotting

generate_data.jl
Transfer data to be csv file, and as input to the R stan model

Nonpooled-Nonstationary-Model
Nonstationary model for individual stations. Stan model and R code to run the model

Pooled-Stationary-Model
Spatially pooled stationary model. Stan model and R code to run the model

Spatially-Varying-Covariate-Model
The final model. Stan model and R code to run the model

validation_metrics.qmd
Different validation metrics for the above three frameworks

validation_plots.qmd
Validation plots including comparison with Atlas 14, cross-validation among different subsets, comparing the three frameworks