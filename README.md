# Average least squares effects

This repo contains reproduction code for simulation and the applied example for the paper "Optimally weighted average derivative effects"
by Oliver Hines, Karla Diaz-Ordaz, and Stijn Vansteelandt

## Simulation study

Simulation code is contained in `R/simulation`. To run simulations, run `experiment_1.R` from the repo root. Simulation data, log files and plots will be saved to `Output/`.

## Illustrated example

Simulation is based on data from the International Warfarin Pharmacogenetics Consortium (IWPC). This data is publicly available under a [Creative Commons Attribution-ShareAlike 4.0 International License](https://www.pharmgkb.org/page/dataUsagePolicy). Note that repo does not contain the IWPC the data, but contains scripts which download it from [PharmaGKB](https://www.pharmgkb.org/downloads).

Analysis code is contained in `R/iwpc`. To load the data, first run `prep_data.R` from the repo root. Then run `analysis.R`. Results are saved to `Output/`.

## Dependencies

This code has been tested using `R v4.2.1`.
Please make sure the following CRAN dependencies are installed:

- `tidyverse`
- `arrow`
- `ranger`
- `mgcv`
- `glmnet`
- `gam`
- `xgboost`
- `SuperLeaner`
- `cowplot`
- `latex2exp`
