# Inverse Probability Weighting-based Mediation Analysis for Microbiome Data
---


## 📖 Overview

This repository contains R scripts for conducting simulation studies under both **null** and **alternative** scenarios to estimate the interventional indirect effect (IIE).
The workflow includes data generation, model fitting, true-value calculation, and bootstrap procedures.

---

## 📂 Repository Structure

| File | Description |
|------|--------------|
| **`null_case.R`** | Runs simulations under null scenarios where indirect effects are absent. |
| **`alternative_case.R`** | Runs simulations under alternative scenarios where indirect effects are present. |
| **`a_simulated_sample_of_X.csv`** | An example for a pseudo-dataset of baseline covariates X, resampled from a real dataset. |
| **`simulate_data.R`** | Generates simulated datasets under specified parameters. |
| **`estimate_true_value_of_IIE.R`** | Computes the true IIE value under alternative scenarios. |
| **`Outcome_model_fit.R`** | Estimates the outcome model. |
| **`estimate_IIE.R`** | Estimates the IIE. |
| **`bootstrap.R`** | Estimates the IIE based on bootstrap samples. |

All code is written in **base R**, using only **CRAN-available packages** for full reproducibility.

---

## 🧠 Scientific Context

This code accompanies the manuscript:

> **Zhang, Y.**, Wang, J., Shen, J., Galloway-Peña, J., Shelburne, S.,  Wang, L., & Hu, J. (2025+).  
> *Inverse Probability Weighting-based Mediation Analysis for Microbiome Data.*  

The framework estimates:
- **Interventional Indirect Effect (IIE)**  

for high-dimensional mediators, while adjust for baseline covariates and an exposure-indecuded mediator-outcome confounder using  **Inverse Probability Weighting** with **SCAD-based regularization**.

---

## ⚙️ System Requirements

### 💻 Software Requirements

| Component | Requirement |
|------------|-------------|
| **R version** | ≥ 4.1 (tested on R 4.1.1) |
| **Operating Systems Tested** | 65-core node equipped with an Intel Cascade Lake CPU |
| **Required Packages** | `snowfall`, `MASS`, `dplyr`, `tidyverse`, `dglm`, `ncvreg`, `glmnet`, `cubature`, `deepTL`|

All packages are platform-independent and available via **CRAN** and **github**.

## 📦 Installation Guide

### 🧰 Package Installation

Before running the R scripts, install required packages from **CRAN** and **github**:

```r
install.packages(c("devtools", "snowfall", "MASS", "dplyr", "tidyverse", "dglm", "ncvreg", "glmnet", "cubature"))

devtools::install_github("SkadiEye/deepTL")
```

  ## 🧭 Instructions for Use

### 1️⃣ Run the Complete Workflow

From your R terminal or RStudio console, execute:

```r
# Step 1. Set working directory to this repository 
setwd("path/to/CausalMediationMicrobiome")

# Step 2. Run one of the main scripts:
source("null_case.R")

source("alternative_case.R")
```
Each script automatically sources all required subfiles and executes the full simulation procedure.

### 2️⃣ Expected Output
Running either main script will:

-  Generate simulated datasets

-  Fit the necessary models

-  Estimate the individual indirect effect (IIE)

-  Save the results to the output directory specified in the script

No additional user steps are needed unless you want to modify default settings.

### 3️⃣ Adjusting Simulation Settings
If you want to explore different simulation configurations (e.g., sample size, number of mediators, signal strength, noise levels),
you can edit the parameter blocks near the top of the main files:
- null_case.R

- alternative_case.R

These files contain clearly labeled parameter sections so you can easily customize the simulation environment.

## 🪪 License

This repository is distributed under the **MIT License**.  
You are free to use, modify, and distribute this code with proper attribution.  
See the [LICENSE](LICENSE) file for full terms.


