# Evaluate a causal estimator against a known DGP

Repeatedly draws datasets from a
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
object, applies a user-supplied estimator to each draw, and returns a
tidy summary of estimator performance against the known ground truth.

## Usage

``` r
causalsim_eval(
  dgp,
  estimator,
  reps = 200L,
  metrics = c("bias", "rmse", "coverage", "power"),
  seed = NULL
)
```

## Arguments

- dgp:

  A `causalsim_dgp` object.

- estimator:

  A function that accepts a data frame produced by
  [`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
  and returns a named numeric vector or single-row data frame. Must
  include an `estimate` field. For `"coverage"` and `"power"` metrics,
  must also include `ci_lower` and `ci_upper`. See Details.

- reps:

  Positive integer. Number of simulation replications. Default `200L`.

- metrics:

  Character vector. Subset of `"bias"`, `"rmse"`, `"coverage"`,
  `"power"`. Default: all four.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the first
  replication. Default `NULL`.

## Value

An S3 object of class `causalsim_eval` with components:

- `metrics`:

  Tidy data frame: `metric`, `value`, `se`.

- `draws`:

  Data frame of per-replication estimates, one row per rep.

- `true_ate`:

  The DGP's true ATE.

- `reps`:

  Number of replications run.

## Details

### Estimator convention

The estimator receives the full data frame returned by
[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md),
including ground-truth columns `.tau` and `.p`. It should return a named
numeric vector or one-row data frame with at minimum:

    my_estimator <- function(data) {
      fit <- lm(Y ~ A + W, data = data)
      est <- coef(fit)["A"]
      se  <- sqrt(vcov(fit)["A", "A"])
      c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
    }

### Metrics

- `"bias"`:

  Monte Carlo estimate of bias: `mean(estimate) - true_ate`.

- `"rmse"`:

  Root mean squared error: `sqrt(mean((estimate - true_ate)^2))`.

- `"coverage"`:

  Proportion of replications where the CI brackets `true_ate`. Requires
  `ci_lower` and `ci_upper`.

- `"power"`:

  Proportion of replications where the CI excludes zero — i.e., the
  rejection rate of H0: ATE = 0. When `true_ate != 0`, this is power;
  when `true_ate == 0`, it is the Type I error rate. Requires `ci_lower`
  and `ci_upper`.

Monte Carlo standard errors (MCSE) are reported alongside each metric.
For bias: `sd(estimates) / sqrt(reps)`. For RMSE: delta-method
approximation. For coverage and power: Bernoulli SE.

### Raw draws

Per-replication estimates are stored in `result$draws` and can be used
for custom analysis or plotting.

## Examples

``` r
dgp <- causalsim_dgp(n = 300, n_confounders = 1, effect = 2,
                     propensity = 0.5)

# OLS estimator (unbiased under RCT)
ols_estimator <- function(data) {
  fit <- lm(Y ~ A + W, data = data)
  est <- coef(fit)[["A"]]
  se  <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

result <- causalsim_eval(dgp, ols_estimator, reps = 100L, seed = 1L)
result
#> <causalsim_eval>  reps: 100  true ATE: 2
#> 
#>    metric   value      se
#>      bias -0.0061 0.01300
#>      rmse  0.1295 0.00723
#>  coverage  0.9300 0.02551
#>     power  1.0000 0.00000
```
