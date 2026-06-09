## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  out.width = "80%",
  warning = FALSE,
  message = FALSE
)
library(rbbnp)
library(ggplot2)
library(gridExtra)

## ----eval=FALSE---------------------------------------------------------------
# biasBound_density(
#   X,                    # Data vector
#   x = NULL,             # Evaluation points (auto if NULL)
#   h = NULL,             # Bandwidth (auto-selected if NULL)
#   h_method = "cv",      # "cv" or "silverman"
#   alpha = 0.05,         # Confidence level
#   kernel.fun = "Schennach2004",
#   xi_lb = NULL,         # Frequency range lower bound
#   xi_ub = NULL          # Frequency range upper bound
# )

## ----basic-density------------------------------------------------------------
# Generate sample data from 2-fold uniform distribution
X <- gen_sample_data(size = 1000, dgp = "2_fold_uniform", seed = 42)

# Estimate density
fit <- biasBound_density(X, h = 0.08, kernel.fun = "Schennach2004")

# Display results
fit

## ----output-structure---------------------------------------------------------
# Key parameters
coef(fit)

# Detailed summary
summary(fit)

## ----density-vis--------------------------------------------------------------
plot(fit)

## ----ft-analysis--------------------------------------------------------------
plot(fit, type = "ft")

## ----cv-bw--------------------------------------------------------------------
fit_cv <- biasBound_density(X, h = NULL, h_method = "cv", kernel.fun = "normal")
cat("CV bandwidth:", coef(fit_cv)["h"], "\n")

plot(fit_cv) + ggtitle("Cross-Validation Bandwidth")

## ----silverman-bw-------------------------------------------------------------
fit_silv <- biasBound_density(X, h = NULL, h_method = "silverman", kernel.fun = "normal")
cat("Silverman bandwidth:", coef(fit_silv)["h"], "\n")

plot(fit_silv) + ggtitle("Silverman's Rule Bandwidth")

## ----direct-bw----------------------------------------------------------------
# Select bandwidth without full estimation
h_cv <- select_bandwidth(X, method = "cv", kernel.fun = "normal")
h_silv <- select_bandwidth(X, method = "silverman", kernel.fun = "normal")

cat("CV:", h_cv, "\nSilverman:", h_silv)

## ----kernel-table, echo=FALSE-------------------------------------------------
knitr::kable(data.frame(
  Kernel = c("Schennach2004", "sinc", "normal", "epanechnikov"),
  Order = c("Infinite", "Infinite", "2", "2"),
  `Best For` = c("High smoothness (recommended)", "Any smoothness level",
                  "Moderate smoothness", "Limited smoothness"),
  check.names = FALSE
))

## ----kernel-comparison, fig.width=9, fig.height=9, out.width="100%"-----------
fit_sch <- biasBound_density(X, kernel.fun = "Schennach2004")
fit_sinc <- biasBound_density(X, kernel.fun = "sinc")
fit_norm <- biasBound_density(X, kernel.fun = "normal")
fit_epan <- biasBound_density(X, kernel.fun = "epanechnikov")

grid.arrange(
  plot(fit_sch) + ggtitle("Schennach2004 (infinite order)"),
  plot(fit_sinc) + ggtitle("Sinc (infinite order)"),
  plot(fit_norm) + ggtitle("Normal (2nd order)"),
  plot(fit_epan) + ggtitle("Epanechnikov (2nd order)"),
  ncol = 2
)

## ----custom-freq, fig.width=6, fig.height=10.5, out.width="100%"--------------
# Default (automatic)
fit_auto <- biasBound_density(X, h = 0.08)

# Custom range
fit_custom <- biasBound_density(X, h = 0.08, xi_lb = 2, xi_ub = 8)

grid.arrange(
  plot(fit_auto, type = "ft") + ggtitle("Automatic Frequency Range"),
  plot(fit_custom, type = "ft") + ggtitle("Custom Range [2, 8]"),
  ncol = 1
)

## ----extract-results----------------------------------------------------------
# Confidence intervals as matrix
ci <- confint(fit)
head(ci, 10)

# Evaluation points
x_points <- fit$x
head(x_points)

# Density estimates
f_hat <- fit$density
head(f_hat)

## ----comparison, fig.width=6, fig.height=10.5, out.width="100%"---------------
h_opt <- select_bandwidth(X, method = "cv", kernel.fun = "sinc")

# Bias-bound with optimal bandwidth
result_opt <- biasBound_density(X, h = h_opt, kernel.fun = "sinc")

# Undersmoothing (half optimal)
result_under <- biasBound_density(X, h = h_opt * 0.5, kernel.fun = "sinc")

grid.arrange(
  plot(result_opt) + ggtitle(paste0("Bias Bound (h = ", round(h_opt, 3), ")")),
  plot(result_under) + ggtitle(paste0("Undersmoothing (h = ", round(h_opt/2, 3), ")")),
  ncol = 1
)

