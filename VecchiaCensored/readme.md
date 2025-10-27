# VecchiaCensored

A comprehensive R package for spatial modeling of censored and non-censored observations using Vecchia approximations, supporting both Gaussian Processes (GP) and Spatially Varying Coefficients (SVC) models in Bayesian and frequentist frameworks.

## Overview

**VecchiaCensored** provides scalable methods for analyzing spatial data with:
- **Censored observations** (e.g., below detection limit measurements)
- **Spatially varying coefficients (SVC)** for modeling non-stationary spatial relationships
- **Vecchia approximation** for computational efficiency with large datasets
- **Bayesian and frequentist estimation** methods
- **Prediction and conditional simulation** capabilities
- **Fixed or varying spatial range parameters** across SVC components

## Key Features

- **Flexible spatial models**: Simple GP, SVC, or hybrid approaches
- **Handles censoring**: Proper treatment of left-censored data
- **Scalable**: Vecchia approximation enables analysis of large spatial datasets
- **Dual frameworks**: Choose between Bayesian (Stan) or frequentist (optimization) estimation
- **Prediction**: Out-of-sample prediction with uncertainty quantification
- **Simulation**: Conditional simulation from fitted models
- **Coefficient prediction**: Predict spatially varying coefficients at new locations

## Installation

### Prerequisites

```r
# Required R packages
install.packages(c("GpGp", "Rfast", "Rcpp", "cmdstanr"))

# For Bayesian estimation, install CmdStan
cmdstanr::install_cmdstan()
```

### Install VecchiaCensored

```r
# From GitHub (once published)
# install.packages("devtools")
devtools::install_github("yacine-idir/VecchiaCensored")
```

### Manual Installation

1. Clone the repository
2. Compile C++ code: `Rcpp::sourceCpp("sparse_start.cpp")`
3. Compile Stan models using `cmdstanr::cmdstan_model()`
4. Source R scripts in order

## Quick Start

### Example 1: Simple GP Model (No Censoring)

```r
library(VecchiaCensored)

# Simulate spatial data
set.seed(123)
n <- 500
locs <- cbind(runif(n), runif(n))
X <- cbind(1, rnorm(n))
Y <- X %*% c(2, 0.5) + rnorm(n)

# Fit model
fit <- fit_model(
  Y_obs = Y,
  locs_obs = locs,
  X_obs = X,
  censored_indices = rep(0, n),  # No censoring
  M = 30,                         # Number of neighbors
  bayesian = FALSE
)

# Predict at new locations
locs_pred <- cbind(runif(50), runif(50))
X_pred <- cbind(1, rnorm(50))

predictions <- prediction(
  fit_object = fit,
  locs_pred = locs_pred,
  X_pred = X_pred,
  M = 30,
  simulations = TRUE,
  n_simulations = 100
)
```

### Example 2: SVC Model with Censoring

```r
# Mark observations below threshold as censored
threshold <- quantile(Y, 0.2)
censored_idx <- ifelse(Y < threshold, 1, 0)
Y_cen <- Y
Y_cen[censored_idx == 1] <- threshold

# Fit SVC model with spatially varying intercept
fit_svc <- fit_model(
  Y_obs = Y_cen,
  locs_obs = locs,
  X_obs = X,
  svc_indices = 1,              # Intercept varies spatially
  censored_indices = censored_idx,
  M = 30,
  bayesian = FALSE,
  fixed_range = TRUE            # Share range parameter across SVCs
)

# Predict including coefficient values
pred_svc <- prediction(
  fit_object = fit_svc,
  locs_pred = locs_pred,
  X_pred = X_pred,
  M = 30,
  pred_coef = TRUE              # Also predict SVC coefficients
)
```

### Example 3: Bayesian Estimation

```r
fit_bayes <- fit_model(
  Y_obs = Y_cen,
  locs_obs = locs,
  X_obs = X,
  svc_indices = c(1, 2),        # Both coefficients vary spatially
  censored_indices = censored_idx,
  M = 30,
  bayesian = TRUE,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  parallel_chains = 4,
  fixed_range = FALSE           # Each SVC has its own range
)

# Access posterior samples
beta_samples <- fit_bayes$parameters$beta
sigma_samples <- fit_bayes$parameters$sigma
```

## Core Functions

### `fit_model()`

Main estimation function supporting multiple model types.

**Arguments:**
- `Y_obs`: Numeric vector of observations
- `locs_obs`: Matrix of spatial coordinates (n × 2)
- `X_obs`: Design matrix of covariates (n × p)
