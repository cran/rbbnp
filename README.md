# rbbnp <img src="man/figures/logo.png" align="right" height="139" alt="rbbnp logo" />

<!-- The R-CMD-check badge is omitted because the development repository is private. -->
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/rbbnp)](https://CRAN.R-project.org/package=rbbnp)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

## Overview

The **rbbnp** package implements the bias-bound approach of [Schennach 2020](https://doi.org/10.1093/restud/rdz065): it bounds the estimation bias via Fourier analysis, giving valid confidence intervals for kernel density and conditional-expectation estimators at optimal, MSE-minimizing bandwidths without undersmoothing.

## Installation

```r
# Install from CRAN
install.packages("rbbnp")

# Or install development version from GitHub
# install.packages("devtools")
devtools::install_github("xinyu-daidai/rbbnp-dev")
```

## Key Functions

| Function | Purpose |
|----------|---------|
| `biasBound_density()` | Density estimation with bias-aware confidence intervals |
| `biasBound_condExpectation()` | Regression with bias-aware confidence intervals |
| `select_bandwidth()` | Cross-validation or Silverman bandwidth selection |

## Usage

### Density Estimation

```r
library(rbbnp)

# Generate sample data
X <- gen_sample_data(size = 500, dgp = "2_fold_uniform", seed = 123)

# Estimate density with bias-aware confidence intervals
fit <- biasBound_density(X, h = 0.1, kernel.fun = "Schennach2004")

# View results
fit
#> Bias-Bound Density Estimation
#> ==============================
#> Observations: 500 | Bandwidth: 0.100 | Kernel: Schennach2004
#> Smoothness: A = 4.30, r = 2.00

# Visualize
plot(fit)
```

### Conditional Expectation (Regression)

```r
# Generate regression data
Y <- -X^2 + 3*X + rnorm(500) * X

# Estimate E[Y|X]
fit_reg <- biasBound_condExpectation(Y, X, h = 0.1)

# Visualize
plot(fit_reg)
```

## Working with Results

Both functions return S3 objects with standard methods:

```r
# Extract parameters (A, r, B, h)
coef(fit)

# Get confidence intervals
confint(fit)

# Detailed summary
summary(fit)

# For regression: fitted values
fitted(fit_reg)
```

## Learning More

- [Get Started](https://xinyudai.net/rbbnp-dev/articles/rbbnp.html) - Quick introduction and basic workflow
- [Density Estimation](https://xinyudai.net/rbbnp-dev/articles/density-estimation.html) - Detailed guide to density estimation
- [Regression](https://xinyudai.net/rbbnp-dev/articles/regression.html) - Conditional expectation estimation
- [Theoretical Background](https://xinyudai.net/rbbnp-dev/articles/theory.html) - Mathematical foundations

## Citation

If you use **rbbnp**, please cite the package (run `citation("rbbnp")` for the current version):

> Dai, X. and Schennach, S. M. (2026). *rbbnp: A Bias Bound Approach to Non-Parametric Inference*. R package version 1.1.0. https://CRAN.R-project.org/package=rbbnp

```bibtex
@Manual{rbbnp,
  title  = {rbbnp: A Bias Bound Approach to Non-Parametric Inference},
  author = {Xinyu Dai and Susanne M. Schennach},
  year   = {2026},
  note   = {R package version 1.1.0},
  url    = {https://CRAN.R-project.org/package=rbbnp},
}
```

The package implements the method introduced in [Schennach (2020)](https://doi.org/10.1093/restud/rdz065).

## Getting Help

- Contact: Xinyu Dai (xinyu_dai@brown.edu)
