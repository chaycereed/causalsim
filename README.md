# causalsim

An R package for defining causal data generating processes with known ground truth and evaluating estimator performance against them.

## Features

- Structural causal model with explicit effect, propensity, and baseline functions
- Named covariate roles: confounder, instrument, effect modifier, noise
- Preset confounding levels (`"low"`, `"moderate"`, `"high"`) or custom functions
- Exact or Monte Carlo true ATE computed at construction time
- Flexible estimator interface: named numeric vector, named list, or one-row data frame
- Tidy performance metrics: bias, RMSE, coverage, and power with Monte Carlo standard errors
- Grid evaluation over the Cartesian product of any DGP parameters
- Reproducible via seed control at every stage

## Installation

Requires R 4.0 or higher. Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("chaycereed/causalsim")
```

## Usage

### Define a data generating process

`causalsim_dgp()` specifies the structural model. Parameters can be scalars, preset strings, or functions of the covariate names.

```r
library(causalsim)

dgp <- causalsim_dgp(
  n             = 500,
  n_confounders = 1,
  effect        = 2,
  propensity    = "moderate",
  baseline      = "moderate"
)
dgp
```

### Draw a dataset

`causalsim_draw()` simulates one dataset from the DGP. The returned data frame includes covariate columns, treatment `A`, outcome `Y`, individual effect `.tau`, and propensity `.p`.

```r
dat <- causalsim_draw(dgp, seed = 1L)
head(dat)
```

### Evaluate an estimator

An estimator is any function that accepts a data frame and returns a named numeric vector with at minimum an `estimate` field. `ci_lower` and `ci_upper` enable coverage and power metrics.

```r
ols_est <- function(data) {
  fit <- lm(Y ~ A + W, data = data)
  est <- coef(fit)[["A"]]
  se  <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

result <- causalsim_eval(dgp, ols_est, reps = 200L, seed = 1L)
result

summary(result)
plot(result)
```

### Evaluate across a parameter grid

`causalsim_grid()` runs the evaluator over the Cartesian product of any DGP parameters, returning a tidy data frame of metrics for each cell.

```r
grid_result <- causalsim_grid(
  dgp       = dgp,
  estimator = ols_est,
  vary      = list(n = c(100L, 250L, 500L, 1000L)),
  reps      = 200L,
  metrics   = c("bias", "rmse"),
  seed      = 1L
)
grid_result
```

## License

MIT License. See `LICENSE` for details.
