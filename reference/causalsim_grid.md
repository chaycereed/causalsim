# Evaluate an estimator across a parameter grid

Runs
[`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
over the Cartesian product of one or more DGP parameter values,
returning a tidy data frame of performance metrics for every cell.
Designed for studying how estimator behavior changes with sample size,
confounding strength, effect size, or noise level.

## Usage

``` r
causalsim_grid(
  dgp,
  estimator,
  vary,
  reps = 200L,
  metrics = c("bias", "rmse", "coverage", "power"),
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- dgp:

  A `causalsim_dgp` object. Provides fixed parameter values for all
  dimensions not listed in `vary`.

- estimator:

  A function as accepted by
  [`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md).

- vary:

  Named list of atomic vectors. Each name must be a valid
  [`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
  argument (except `covariates`, which cannot be varied atomically).
  Each vector supplies the values to try for that dimension. The full
  grid is the Cartesian product of all dimensions.

- reps:

  Positive integer. Replications per grid cell. Default `200L`.

- metrics:

  Character vector. Passed to
  [`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md).

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the first
  cell. The same RNG stream continues across cells, so results are
  jointly reproducible. Default `NULL`.

- verbose:

  Logical. If `TRUE`, prints a progress message before each cell.
  Default `FALSE`.

## Value

An S3 object of class `causalsim_grid` with components:

- `results`:

  Tidy data frame: one row per (cell, metric). Columns: grid parameter
  names, `metric`, `value`, `se`.

- `grid`:

  Data frame of grid cells, one row per cell.

- `vary`:

  Character vector of varied parameter names.

- `reps`:

  Replications per cell.

- `metrics`:

  Metrics evaluated.

## Details

### How it works

For each cell in `expand.grid(vary)`, `causalsim_grid()` takes `dgp`'s
stored original parameters, overrides the cell's values, reconstructs a
new
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md),
runs
[`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md),
and tags the resulting metrics with the cell's parameter values.

### Constraints on `vary`

Elements of `vary` must be atomic vectors (character, numeric, integer).
Functions cannot be varied via this interface. `covariates` (a list of
[`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
objects) is excluded. For complex covariate variation, construct DGPs
manually and use
[`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
directly.

Note that varying `n_confounders` from 1 to 2 changes auto-generated
covariate names from `W` to `W1, W2`. If the base DGP's `effect`
function references `W`, the reconstructed DGP will error at
construction time. Design the base DGP accordingly.

## Examples

``` r
dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2,
                     propensity = 0.5)

ols_estimator <- function(data) {
  fit <- lm(Y ~ A + W, data = data)
  est <- coef(fit)[["A"]]
  se  <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

# \donttest{
grid_result <- causalsim_grid(
  dgp = dgp,
  estimator = ols_estimator,
  vary = list(n = c(250L, 500L, 1000L)),
  reps = 50L,
  seed = 1L
)
grid_result
#> <causalsim_grid>  3 cells  vary: n  reps/cell: 50
#>   metrics: bias, rmse, coverage, power
#> 
#>     n   metric        value          se
#>   250     bias -0.018135254 0.019081996
#>   250     rmse  0.134799456 0.011981623
#>   250 coverage  0.900000000 0.042426407
#>   250    power  1.000000000 0.000000000
#>   500     bias -0.007350673 0.012609789
#>   500     rmse  0.088574059 0.009552722
#>   500 coverage  0.940000000 0.033585711
#>   500    power  1.000000000 0.000000000
#>  1000     bias  0.010538141 0.010225101
#>  1000     rmse  0.072347315 0.006066306
#> 
#>   ... 2 more rows. Access full results via $results.
# }
```
