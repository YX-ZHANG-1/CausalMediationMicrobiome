# Inverse Probability Weighting-Based Mediation Analysis for Microbiome Data

## Overview

This repository contains the R code used for simulation studies of an inverse probability weighting (IPW) estimator of the interventional indirect effect (IIE) with high-dimensional compositional mediators. The simulations include a binary exposure, a binary exposure-induced mediator-outcome confounder, centered log-ratio (clr) transformed mediators, two baseline covariates, and a binary outcome.

Two simulation scenarios are provided:

- `null_case.R` evaluates estimation performance and type I error when the true IIE is zero.
- `alternative_case.R` evaluates estimation performance and power when the true IIE is nonzero. It first approximates the true IIE using large simulated samples.

The default implementation uses SCAD-penalized outcome regression with AIC tuning, a deep-learning ensemble for the exposure-induced mediator-outcome confounder model, IPW estimation, normal-approximation confidence intervals based on bootstrap standard errors, and parallel computation.

## Repository contents

| File | Description |
| --- | --- |
| [`null_case.R`](null_case.R) | Main script for the null scenario. It reports Monte Carlo bias, standard deviation, RMSE, and type I error. |
| [`alternative_case.R`](alternative_case.R) | Main script for the alternative scenario. It approximates the true IIE and reports Monte Carlo bias, standard deviation, RMSE, and power. |
| [`simulate_data.R`](simulate_data.R) | Generates `Y`, `Z`, `L`, two baseline covariates, and clr-transformed compositional mediators. |
| [`Outcome_model_fit.R`](Outcome_model_fit.R) | Fits penalized outcome models and selects the tuning parameter by cross-validation, AIC, BIC, or extended BIC. |
| [`estimate_IIE.R`](estimate_IIE.R) | Estimates the IIE using the proposed IPW procedure. |
| [`bootstrap.R`](bootstrap.R) | Generates bootstrap IIE estimates. |
| [`estimate_true_value_of_IIE.R`](estimate_true_value_of_IIE.R) | Approximates the true IIE under the alternative scenario using large simulated samples and numerical integration. |
| [`a pseudo-dataset for baseline covariates X.csv`](a%20pseudo-dataset%20for%20baseline%20covariates%20X.csv) | Pseudo-data containing `age` and binary `sex`, used in place of the real baseline covariate data. |

## Simulation workflow

Each Monte Carlo replication performs the following steps:

1. Resample the two baseline covariates from the supplied pseudo-dataset.
2. Simulate the binary exposure `Z`.
3. Simulate the binary exposure-induced mediator-outcome confounder `L`.
4. Simulate the clr-transformed compositional mediators `clrM`.
5. Simulate the binary outcome `Y`.
6. Fit a penalized outcome model. The exposure, confounder, and baseline covariates are not penalized; the mediator coefficients are penalized.
7. Estimate the conditional models required for the IPW estimator.
8. Estimate the IIE and obtain bootstrap estimates.
9. Summarize Monte Carlo bias, standard deviation, RMSE, and either type I error or power.

The active scripts assume exactly two baseline covariates. In the supplied pseudo-dataset, the first covariate is age and is divided by 100 before analysis; the second is binary sex.

## Software requirements

The code requires R and the following packages:

- `snowfall`
- `MASS`
- `dplyr`
- `tidyr`
- `tibble`
- `dglm`
- `ncvreg`
- `glmnet`
- `cubature`
- `deepTL`

The `cubature` package is required by the alternative-scenario true-IIE calculation. The `deepTL` package supplies `importDnnet()`, `ensemble_dnnet()`, and the associated prediction method used by the IPW estimator.

Install the CRAN packages with:

```r
install.packages(
  c(
    "snowfall",
    "MASS",
    "dplyr",
    "tidyr",
    "tibble",
    "dglm",
    "ncvreg",
    "glmnet",
    "cubature",
    "devtools"
  )
)
```

Install `deepTL` from GitHub with:

```r
devtools::install_github("SkadiEye/deepTL")
```

## Running the simulations

Clone or download the repository, open a terminal in the repository directory, and run one of the main scripts.

Null scenario:

```sh
Rscript null_case.R
```

Alternative scenario:

```sh
Rscript alternative_case.R
```

Alternatively, from R or RStudio:

```r
setwd("/path/to/CausalMediationMicrobiome")

source("null_case.R")
# or
source("alternative_case.R")
```

The working directory must contain the main script, the supporting R files, and the pseudo-dataset because the scripts use relative paths. Run the null and alternative scripts in separate R sessions because each script begins by clearing the current workspace.

## Default settings

| Setting | Null scenario | Alternative scenario |
| --- | ---: | ---: |
| Monte Carlo replications (`nrep`) | 500 | 500 |
| Bootstrap replications (`nbts`) | 400 | 400 |
| Subjects per simulated dataset (`n_subj`) | 70 | 70 |
| Mediator components (`n_M`) | 134 | 134 |
| Large-sample size for true-IIE approximation (`n_true`) | Not applicable | 10,000 |
| Maximum attempts per estimate (`iter1`) | 50 | 50 |
| Maximum accepted absolute IIE (`threshold`) | 1 | 1 |
| Penalization | SCAD | SCAD |
| Tuning method | AIC | AIC |
| Parallel workers (`ncpus`) | 40 | 40 |

These defaults are computationally intensive. For a short installation test, temporarily reduce `nrep`, `nbts`, and `ncpus`. In the alternative-scenario script, also reduce `n_true`. Do not set `ncpus` above the number of CPU cores allocated to the R job.

## Output

Each main script saves an `.RData` file in the working directory. The filename follows this pattern:

```text
results_penalty_<penalty>_tune_<tune>_case<case>_n<n_subj>_nrep<nrep>.RData
```

For example:

```text
results_penalty_SCAD_tune_aic_casenull_n70_nrep500.RData
results_penalty_SCAD_tune_aic_casealternative_n70_nrep500.RData
```

The null-scenario output contains:

- `true.parameters`: simulation parameters
- `IIE.true`: true IIE, fixed at zero
- `IIE_hat_all`: Monte Carlo IIE estimates
- `IIE_hat_bs_all`: bootstrap estimates for each Monte Carlo sample
- `IIE_bias`, `IIE_sd`, and `IIE_rmse`: Monte Carlo performance summaries
- `type_one_error`: estimated type I error based on normal bootstrap intervals
- `elapsed_time`: elapsed runtime in seconds

The alternative-scenario output contains analogous estimates and performance summaries, with `power` in place of `type_one_error`. It also contains:

- `IIE_true_all`: large-sample estimates used to approximate the true IIE
- `power`: estimated rejection probability for the null hypothesis that the IIE is zero

Load a saved result with:

```r
load("results_penalty_SCAD_tune_aic_casenull_n70_nrep500.RData")
```

## Customizing the simulation

The main scripts contain labeled parameter blocks near the top of each file. Common settings to modify include:

- `nrep`: number of Monte Carlo replications
- `nbts`: number of bootstrap replications
- `n_subj`: sample size
- `n_M`: number of mediator components
- `n_true`: large-sample size for approximating the true IIE in the alternative scenario
- `ncpus`: number of parallel workers
- `penalty`: `"SCAD"`, `"MCP"`, or `"lasso"`
- `tune`: `"cv"`, `"aic"`, `"bic"`, or `"ebic"`
- `iter1`: maximum number of attempts after estimation failure
- `threshold`: maximum accepted absolute IIE estimate
- regression coefficients governing the exposure, confounder, mediator, and outcome models

If the baseline-covariate file is replaced, retain two numeric columns or update the column-selection and model code accordingly. Ensure that `n_M` equals the length of `beta_M`.

## Reproducibility notes

- The simulation and resampling steps use explicitly constructed random-number seeds.
- The scripts save the parameter settings and elapsed runtime with the results.
- Exact numerical results can still depend on the R version, package versions, operating system, and parallel-computing environment.
- Failed model fits and inadmissible IIE estimates are retried up to `iter1` times. If no valid estimate is obtained, the corresponding result is stored as `NA`.
- The alternative scenario is substantially more computationally intensive because it uses large samples and numerical integration to approximate the true IIE.

## Citation

If you use this code, please cite the accompanying manuscript:

> Zhang, Y., Wang, J., Shen, J., Galloway-Peña, J., Shelburne, S., Wang, L., and Hu, J. *Inverse Probability Weighting-Based Mediation Analysis for Microbiome Data.* Manuscript submitted to *Bioinformatics*.

After the software is archived on Zenodo, add the version-specific software DOI and citation here.

## License

This software is distributed under the [MIT License](LICENSE.txt). Copyright (c) 2025 Yuexia Zhang.
