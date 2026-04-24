# A Complete Simulation Study with causalsim

## Overview

Evaluating causal estimators requires knowing the answer. `causalsim`
gives you that by making the data-generating process explicit: you
specify the structural model, the package simulates data from it, and
you measure how well any estimator recovers the truth you specified.

This vignette walks through a complete simulation study:

1.  Define a DGP with one confounder and a known treatment effect
2.  Inspect a draw from it to confirm the structural model
3.  Compare a naive (unadjusted) estimator to OLS adjustment
4.  Vary sample size with
    [`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
    to study finite-sample precision
5.  Vary confounding strength to see where naive estimation breaks down

``` r
library(causalsim)
```

## Step 1: Define the DGP

The structural model for this study is:

$$W \sim N(0,1),\quad A \mid W \sim \text{Bernoulli}\!\left( \text{logistic}(0.5\, W) \right),\quad Y = 2A + 0.5W + \varepsilon,\quad\varepsilon \sim N(0,1)$$

In
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
terms: one standard-normal confounder, a constant effect of 2, moderate
propensity confounding, and a moderate baseline shift.

``` r
dgp <- causalsim_dgp(
  n             = 500,
  n_confounders = 1,
  effect        = 2,
  propensity    = "moderate",
  baseline      = "moderate"
)
dgp
#> <causalsim_dgp>
#>   n            : 500
#>   true ATE     : 2.0000
#>   heterogeneous: FALSE
#>   sigma        : 1.00
#>   covariates   :
#>     W       normal  [confounder]
```

The true ATE is exact because `effect` is a scalar — no Monte Carlo
approximation needed.

## Step 2: Inspect a Draw

[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
simulates one dataset. The columns `.tau` and `.p` are the individual
causal effect and propensity score — diagnostic metadata not available
in real observational data.

``` r
dat <- causalsim_draw(dgp, seed = 1L)
head(dat)
#>            W A          Y .tau        .p
#> 1 -0.6264538 0 -1.7679183    2 0.4223273
#> 2  0.1836433 0 -0.7538327    2 0.5229393
#> 3 -0.8356286 0 -1.6682940    2 0.3970399
#> 4  1.5952808 0  1.4649285    2 0.6894695
#> 5  0.3295078 1  0.8739842    2 0.5410956
#> 6 -0.8204684 0 -2.4452377    2 0.3988560
```

With moderate confounding, the treated and control groups differ on the
pre-treatment covariate:

``` r
aggregate(W ~ A, data = dat, FUN = mean)
#>   A          W
#> 1 0 -0.2440553
#> 2 1  0.3162376
```

That difference is exactly what makes naive regression biased — any
estimator that ignores `W` will confuse baseline variation with
treatment effect.

## Step 3: Define Estimators

Both estimators return a named numeric vector with `estimate`,
`ci_lower`, and `ci_upper`. The `ci_*` fields enable coverage and power
metrics in addition to bias and RMSE.

``` r
# Naive: regresses Y on A only — omits the confounder
naive_est <- function(data) {
  fit <- lm(Y ~ A, data = data)
  est <- coef(fit)[["A"]]
  se  <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}

# OLS: adjusts for the observed confounder W
ols_est <- function(data) {
  fit <- lm(Y ~ A + W, data = data)
  est <- coef(fit)[["A"]]
  se  <- sqrt(vcov(fit)["A", "A"])
  c(estimate = est, ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se)
}
```

## Step 4: Evaluate Each Estimator

[`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
runs `reps` independent replications and returns bias, RMSE, coverage,
and power with Monte Carlo standard errors.

``` r
eval_naive <- causalsim_eval(dgp, naive_est, reps = 300L, seed = 1L)
eval_naive
#> <causalsim_eval>  reps: 300  true ATE: 2
#> 
#>    metric  value      se
#>      bias 0.2353 0.00587
#>      rmse 0.2562 0.00571
#>  coverage 0.3700 0.02787
#>     power 1.0000 0.00000
```

``` r
eval_ols <- causalsim_eval(dgp, ols_est, reps = 300L, seed = 1L)
eval_ols
#> <causalsim_eval>  reps: 300  true ATE: 2
#> 
#>    metric   value      se
#>      bias -0.0002 0.00529
#>      rmse  0.0915 0.00338
#>  coverage  0.9433 0.01335
#>     power  1.0000 0.00000
```

The naive estimator’s bias is substantial: treatment is positively
correlated with $W$, which also raises $Y$ through the baseline, so the
unadjusted coefficient absorbs part of that association. OLS eliminates
the bias by conditioning on $W$. Coverage for the naive estimator is
well below the nominal 95 % because confidence intervals are centered on
the wrong value.

## Step 5: Vary Sample Size with `causalsim_grid()`

[`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
evaluates an estimator over the Cartesian product of the supplied
parameter values. Here we vary `n` across four levels to track how the
OLS estimator’s precision improves with more data.

``` r
grid_n <- causalsim_grid(
  dgp       = dgp,
  estimator = ols_est,
  vary      = list(n = c(100L, 250L, 500L, 1000L)),
  reps      = 300L,
  metrics   = c("bias", "rmse"),
  seed      = 1L
)
grid_n
#> <causalsim_grid>  4 cells  vary: n  reps/cell: 300
#>   metrics: bias, rmse
#> 
#>     n metric        value          se
#>   100   bias -0.018028707 0.012266598
#>   100   rmse  0.212874118 0.009119602
#>   250   bias  0.014982725 0.007414731
#>   250   rmse  0.129085151 0.004892668
#>   500   bias  0.003375758 0.005195662
#>   500   rmse  0.089904798 0.003940619
#>  1000   bias  0.003163125 0.003734224
#>  1000   rmse  0.064648194 0.002489179
```

RMSE roughly halves as $n$ quadruples — consistent with $\sqrt{n}$-rate
convergence for OLS in a correctly specified model. Bias stays near zero
at every sample size.

## Step 6: Vary Confounding Strength

We can vary the propensity string to study how bias grows as confounding
increases. Because
[`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
accepts one estimator at a time, we run it separately for each and
combine the results.

``` r
conf_levels <- list(propensity = c("low", "moderate", "high"))

grid_naive <- causalsim_grid(dgp, naive_est,
                              vary    = conf_levels,
                              reps    = 300L,
                              metrics = "bias",
                              seed    = 1L)

grid_ols   <- causalsim_grid(dgp, ols_est,
                              vary    = conf_levels,
                              reps    = 300L,
                              metrics = "bias",
                              seed    = 1L)

comparison <- rbind(
  cbind(estimator = "naive", grid_naive$results),
  cbind(estimator = "ols",   grid_ols$results)
)
comparison <- comparison[order(comparison$propensity, comparison$estimator), ]
rownames(comparison) <- NULL
comparison
#>   estimator propensity metric         value          se
#> 1     naive       high   bias  0.4089432973 0.006098827
#> 2       ols       high   bias -0.0001855971 0.005888572
#> 3     naive        low   bias  0.1305019380 0.005990692
#> 4       ols        low   bias  0.0045047637 0.005556133
#> 5     naive   moderate   bias  0.2419929827 0.005443921
#> 6       ols   moderate   bias -0.0008346139 0.005360097
```

Naive bias grows proportionally with confounding strength. OLS remains
near zero across all three levels because $W$ is observed and included
in the model.

------------------------------------------------------------------------

This three-step workflow — define, evaluate, grid — scales naturally to
more complex settings: heterogeneous effects, multiple confounders, or
estimators that require a different adjustment strategy. See
[`?causalsim_dgp`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
for the full covariate specification API and
[`?covar`](https://chaycereed.github.io/causalsim/reference/covar.md)
for defining non-normal covariates.
