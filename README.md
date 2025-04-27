# Causal Impact of Zero Lower Bound Monetary Policy
## Overview

This repository contains codes and scripts to evaluate the causal effects of ZLB at the state level in the United States. The experimental framework combines BVAR, SCM and ML with block bootstrap inference to construct baseline moderator (G) to split the sample.

## Dependencies
1.**R (>= 4.0)** with packages
  - `vars `, `BVAR`, `Synth`, `boot`, `tidyverse`

2.**Python (>= 3.8)** with packages
- `pandas`, `numpy`, `scikit-learn`, `xgboost`, `matplotlib`
3. **Jupyter** for running notebooks.

## Workflow
1. Data Preprocessing (`preprocessing.ipynb`, `new_preprocessing.ipynb`, `RData_Cleaning.R`)
   - Import state-level time series for GDP, shadow rate, financial indicators.
   - Clean, transform, and merge into panel data format.
2. Bayesian VAR & Feature Extraction (`bayesianVAR.R`)
   - Fit BVAR(p) model in the pre-treatment time period for each state.
   - Extract vectorized VAR coefficients and impulse-response functions (IRFs).
3. Construct Composite Score (G) (`composite_score_supervised.ipynb`, `G-Selection_supervised.ipynb`)
   - Train supervised ML models (Decision Tree, Gaussian Process, Random Forest proximity, XGBoost) to predict pre-treatment GDP growth.
   - Combine predicted scores into a composite score.
   - Choose threshold to split states into group G=1 vs G=0.
4. Synthetic Control & Bootstrap (`synthetic_control_supervised.R`, `synthetic_control_plus_bootstrap.R`)
   - For each treated state, construct synthetic control using SCM with combined dynamic features.
   - Compute treatment effects and block-bootstrap confidence intervals (500 replicates, block length 2).
5. Analysis & Interpretation
   - Aggregate point estimates and confidence intervals across states.
   - Evaluate average treatment effect and heterogeneous effects by G-group.

## Usage
1. **Notebooks:** Launch Jupyter and open `.ipynb` files in sequence.
2. **R Scripts:** Run in RStudio.



