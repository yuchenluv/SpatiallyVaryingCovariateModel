# Spatially Varying Covariate Model
This repository includes the codes for the manuscript [Hierarchical Bayesian Spatial Framework Reveals Robust Increasing Trends in Daily Extreme Rainfall Across the Western Gulf Coast](https://arxiv.org/abs/2502.02000). The paper is currently under review, so this repository is preliminary and subject to updates.

## How to cite

To cite our results/methods, please cite it like follows:

```
@misc{lu2025bayesianspatiotemporalnonstationarymodel,
      title={Bayesian Spatiotemporal Nonstationary Model Quantifies Robust Increases in Daily Extreme Rainfall Across the Western Gulf Coast}, 
      author={Yuchen Lu and Ben Seiyon Lee and James Doss-Gollin},
      year={2025},
      eprint={2502.02000},
      archivePrefix={arXiv},
      primaryClass={stat.AP},
      url={https://arxiv.org/abs/2502.02000}, 
}
```

## How to use

This repository contains all the code and data used in our paper, based entirely on publicly available datasets.

To reproduce our analysis, follow the steps below:

1. **Data Preprocessing**  
   Run `generate_data.jl` to preprocess annual maximum rainfall records and the CO₂ covariate. This will generate `.csv` files that serve as inputs to the statistical models.

2. **Model Estimation**  
   Choose one of the following modeling frameworks:
   - **Pooled Stationary Model**
   - **Nonpooled Nonstationary Model**
   - **Spatially Varying Covariate Model**  
   
   Use the corresponding `.qmd` file to estimate distribution parameters. The estimated parameters will be saved as `.csv` files.

3. **Result Visualization**  
   Use `plot_results.qmd` to generate most of the plots used in the paper.

For more details, please refer to the descriptions of individual files in this repository.

## File Descriptions

All core scripts are located in the `scripts/` directory. Below is a brief description of each file:

### Data Processing
- `get_data.jl`  
  Functions to process annual maximum series (AMS) precipitation from the GHCN dataset and extract climate covariates (e.g., CO₂).

- `params.jl`  
  Defines constraints for extracting AMS precipitation.

- `get_Atlas14.jl`  
  Web scraper to collect Atlas 14 IDF estimates for comparison.

- `generate_data.jl`  
  Converts processed data into `.csv` files, which serve as inputs to the R/Stan models.

### Utilities
- `cal_util.jl`  
  Utility functions for transforming estimated distribution parameters into return level estimates, and computing validation metrics.

- `plot_util.jl`  
  Utility functions for generating plots.

### Modeling Frameworks
- `Nonpooled-Nonstationary-Model/`  
  Implements a nonstationary model for individual stations using Stan and R.

- `Pooled-Stationary-Model/`  
  Implements a spatially pooled stationary model using Stan and R.

- `Spatially-Varying-Covariate-Model/`  
  Implements the final model, with spatially varying covariate effects. Includes Stan and R code.

### Quarto Documents
- `plot_results.qmd`  
  Processes simulation outputs and generates the figures shown in the paper.

- `validation_metrics.qmd`  
  Computes different validation metrics across the three modeling frameworks.

- `validation_plots.qmd`  
  Creates validation plots, including:
  - Comparisons with Atlas 14
  - Cross-validation among data subsets
  - Comparison across modeling frameworks
