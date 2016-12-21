# Time-dependent gene expression Analysis (TDGEA)

These scripts analyze Affymetrix microarray data to look at relationships between time and expression levels.

## Running the Analysis

Generally follows this procedure:
1. Open R
1. Set working directory as `src`
1. Load data with `source("data/load_meta.R");`
1. Create regression models with `source("analysis/probe_regression.R");`
1. Create plots with `source("visualization/plots.R");`