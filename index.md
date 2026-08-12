# causalsim

An R package for defining causal data generating processes with known
ground truth and evaluating estimator performance against them.

## Features

- Structural causal model with explicit effect, propensity, and baseline
  functions
- Named covariate roles: confounder, effect modifier, noise
- Preset confounding levels (`"low"`, `"moderate"`, `"high"`) or custom
  functions
- Exact or Monte Carlo true ATE computed at construction time
- Flexible estimator interface: named numeric vector, named list, or
  one-row data frame
- Tidy performance metrics: bias, RMSE, coverage, and power with Monte
  Carlo standard errors
- Grid evaluation over the Cartesian product of any DGP parameters
- Reproducible via seed control at every stage

## Installation

Requires R 4.0 or higher. Install from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("chaycereed/causalsim")
```

## Usage

### Quick start

[`causalsim()`](https://chaycereed.github.io/causalsim/reference/causalsim.md)
generates a dataset in one call. The returned data frame includes
covariate columns, treatment `A`, outcome `Y`, individual effect `.tau`,
and true propensity `.p`.

``` r

library(causalsim)

data <- causalsim(n = 500, n_confounders = 1, effect = 2, seed = 1L)
head(data)
```

For estimator benchmarking or multi-draw workflows, use
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
and
[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
directly (see below).

### Define a data generating process

[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
specifies the structural model. Parameters can be scalars, preset
strings, or functions of the covariate names.

``` r

dgp <- causalsim_dgp(
  n = 500,
  n_confounders = 1,
  effect = 2,
  propensity = "moderate",
  baseline = "moderate"
)
dgp
```

### Draw a dataset

[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
simulates one dataset from the DGP. The returned data frame includes
covariate columns, treatment `A`, outcome `Y`, individual effect `.tau`,
and propensity `.p`.

``` r

dat <- causalsim_draw(dgp, seed = 1L)
head(dat)
```

### Evaluate an estimator

An estimator is any function that accepts a data frame and returns a
named numeric vector with at minimum an `estimate` field. `ci_lower` and
`ci_upper` enable coverage and power metrics.

``` r

ols_est <- function(data) {
  fit <- lm(Y ~ A + W, data = data)
  est <- coef(fit)[["A"]]
  se <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

result <- causalsim_eval(dgp, ols_est, reps = 200L, seed = 1L)
result

summary(result)
plot(result)
```

### Evaluate across a parameter grid

[`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
runs the evaluator over the Cartesian product of any DGP parameters,
returning a tidy data frame of metrics for each cell.

``` r

grid_result <- causalsim_grid(
  dgp = dgp,
  estimator = ols_est,
  vary = list(n = c(100L, 250L, 500L, 1000L)),
  reps = 200L,
  metrics = c("bias", "rmse"),
  seed = 1L
)
grid_result
```

### Heterogeneous effects

Declare a covariate with `role = "effect_modifier"` and reference it in
a function passed to `effect` to make the treatment effect vary across
subgroups. The individual effects are stored in the `.tau` column as
ground truth.

``` r

het_dgp <- causalsim_dgp(
  n = 2000,
  covariates = list(
    W = covar("normal", role = "confounder"),
    V = covar("binary", role = "effect_modifier", prob = 0.5)
  ),
  effect     = function(V) 2 + 3 * V   # effect is 2 when V = 0, 5 when V = 1
)

d <- causalsim_draw(het_dgp, seed = 1L)
tapply(d$.tau, d$V, mean)              # 2 for V = 0, 5 for V = 1
```

The function passed to `effect` is what activates the modifier; a
covariate labelled `effect_modifier` that `effect` never references is
inert, and
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
warns when that happens.

## Roadmap

Planned for a future release:

- **Additional covariate roles.** The current roles (`confounder`,
  `effect_modifier`, `noise`) cover the common cases. Two single-path
  roles would complete the treatment/outcome taxonomy:
  - `instrument` — drives treatment only (enters the propensity model,
    excluded from the outcome), for benchmarking instrumental-variable
    estimators.
  - `prognostic` — drives the outcome only (enters the baseline,
    independent of treatment), for studying precision covariates and
    variance reduction.

  These will ship together with worked examples that demonstrate each
  (an IV estimator and a variance-reduction comparison, respectively).
- **Assumption-violation helpers.** First-class support for unmeasured
  confounding and positivity violations.

## License

MIT License. See `LICENSE` for details.
