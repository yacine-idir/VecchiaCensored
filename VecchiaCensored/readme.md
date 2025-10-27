# VecchiaCensored

Spatial modeling for censored and non-censored data using Vecchia approximations.

## Features

- Handles censored observations
- Spatially Varying Coefficients (SVC)
- Bayesian and frequentist estimation
- Scalable via Vecchia approximation

## Installation

```r
# Install from GitHub
devtools::install_github("yacine-idir/VecchiaCensored")

# Install CmdStan for Bayesian models
cmdstanr::install_cmdstan()
```

## Quick Example

```r
# Fit a simple GP model
fit <- fit_model(
  Y_obs = Y,
  locs_obs = locs,
  X_obs = X,
  censored_indices = rep(0, length(Y)),
  M = 30,
  bayesian = FALSE
)

# Make predictions
pred <- prediction(fit, locs_pred, X_pred, M = 30)
```

## Documentation

Detailed documentation coming soon.

## License

MIT
