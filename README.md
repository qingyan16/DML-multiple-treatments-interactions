# DML-multiple-treatments-interactions

Code for simulation studies and an HIV-kidney disease application related to double machine learning with multiple treatments and treatment interactions.

## Repository Structure

- `0_data_generation.R`: simulation data generation code.
- `1_PLR_functions.R`: helper functions for the partially linear regression (PLR) model, used by simulations and the application analysis.
- `2_PLR_simulation.R`: simulation script for PLR models with multiple concurrent treatments.
- `3_IRM_simulation.R`: simulation script for interactive regression / multi-valued treatment regimens.
- `4_HIV_application.R`: application script for the HIV-kidney disease analysis.

## Main Requirements

These scripts rely on a mix of R packages for simulation, estimation, and machine learning. The package sets differ slightly by script, but the codebase uses the following packages:

```r
foreach
doParallel
openxlsx
doSNOW
progress
clusterGeneration
mvtnorm
caret
tidyverse
purrr
PSweight
reshape2
gtsummary
mice
glmnet
randomForest
xgboost
ranger
gbm
nnet
DoubleML
mlr3
reticulate
```

`3_IRM_simulation.R` also calls Python modules through `reticulate`, including `sklearn` and `doubleml`.

## How the Code Is Organized

Typical workflow:

1. Run `0_data_generation.R` to generate simulation data.
2. Use `1_PLR_functions.R` as the helper function file for PLR estimation.
3. Run `2_PLR_simulation.R` for PLR simulation experiments.
4. Run `3_IRM_simulation.R` for simulations with multi-valued treatment regimens.
5. Run `4_HIV_application.R` for the HIV-kidney disease application.

## Notes Before Running

- Some scripts use hard-coded local paths through `dir_path`; update those paths to match your environment before running.
- `3_IRM_simulation.R` includes a hard-coded Python virtual environment path in `reticulate`; that path must also be updated locally if you want to reproduce the IRM analysis.
- The simulation code is designed for multi-core or cluster execution and may take substantial time to complete.

## Outputs

The simulation scripts save intermediate or final `.rds` objects to result folders referenced by `dir_path`. Make sure those directories exist before execution.

## Reproducibility

- Seeds are set inside the scripts for replication.
- Because several scripts depend on external file paths and local environments, exact reproduction requires aligning the directory structure and Python/R package setup with the original analysis environment.
