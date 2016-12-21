# Time-dependent gene expression Analysis (TDGEA)

These scripts analyze Affymetrix microarray data to look at relationships between time and expression levels.

## Running the Analysis

Generally follows this procedure:
1. Open R
1. Set working directory as `src`
1. Load data with `source("data/load_meta.R");`
1. Analyze data using a script of choice in the `analysis` folder
1. Load plotting functions with `source("visualization/plotting_functions.R");`
1. Create plots with a corresponding script in the `visualization` folder