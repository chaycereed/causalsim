# Simulate a causal dataset with known ground truth

A convenience wrapper around
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
and
[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md).
Defines a causal data generating process and returns one simulated
dataset in a single call.

## Usage

``` r
causalsim(
  n,
  effect = 1,
  propensity = "moderate",
  baseline = 0,
  sigma = 1,
  covariates = list(),
  n_confounders = 0L,
  n_effect_modifiers = 0L,
  n_instruments = 0L,
  n_noise = 0L,
  mc_draws = 10000L,
  seed = NULL
)
```

## Arguments

- n:

  Positive integer. Sample size for each simulated dataset.

- effect:

  Numeric scalar or function. A scalar specifies a constant
  (homogeneous) treatment effect; ATE = CATE everywhere. A function
  should accept named arguments matching covariate names defined in the
  DGP and return a numeric vector of individual-level causal effects
  (CATE). See Details.

- propensity:

  Numeric scalar, preset string, or function. Treatment assignment
  probability. A scalar (e.g. `0.5`) gives a constant propensity
  (randomized trial). Preset strings `"low"`, `"moderate"`, `"high"`
  generate a logistic propensity over confounders with coefficients
  0.25, 0.5, and 1.0 respectively. A function follows the same
  named-argument convention as `effect` and must return values in \[0,
  1\]. Defaults to `"moderate"`.

- baseline:

  Numeric scalar, preset string, or function. Mean potential outcome
  under control, `E[Y(0) | W]`. Preset strings follow the same levels as
  `propensity` and apply a linear combination of confounders. Defaults
  to `0`.

- sigma:

  Positive numeric. Standard deviation of the outcome noise term.
  Default `1`.

- covariates:

  Named list of
  [`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
  objects (Option A / explicit path). Each name becomes the column name
  in generated data and the argument name expected by `effect`,
  `propensity`, and `baseline` functions. Merged with any auto-generated
  covariates; name collisions error.

- n_confounders:

  Non-negative integer. Standard normal confounders auto-generated as
  `W` (single) or `W1, W2, ...` (multiple).

- n_effect_modifiers:

  Non-negative integer. Auto-generates standard normal effect modifiers
  as `V` or `V1, V2, ...`.

- n_instruments:

  Non-negative integer. Auto-generates standard normal instruments as
  `Z` or `Z1, Z2, ...`.

- n_noise:

  Non-negative integer. Auto-generates standard normal noise covariates
  as `X` or `X1, X2, ...`.

- mc_draws:

  Positive integer. Monte Carlo draws for true ATE approximation.
  Default `10000L`. Ignored for scalar `effect`.

- seed:

  Integer or `NULL`. Passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing.
  Default `NULL`.

## Value

A data frame with covariate columns, `A`, `Y`, `.tau`, and `.p`.

## Details

The returned data frame includes all covariate columns, the treatment
indicator `A`, the outcome `Y`, and two ground-truth columns:

- `.tau`:

  Individual-level causal effect (CATE).

- `.p`:

  True propensity score.

For full control over the DGP (multiple draws, grid evaluation,
estimator benchmarking), use
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
and
[`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
directly.

## Examples

``` r
# One confounder, constant effect of 2, moderate confounding
data <- causalsim(n = 500, n_confounders = 1, effect = 2, seed = 1L)
head(data)
#>            W A          Y .tau        .p
#> 1 -0.6264538 0 -1.4546914    2 0.4223273
#> 2  0.1836433 0 -0.8456543    2 0.5229393
#> 3 -0.8356286 0 -1.2504797    2 0.3970399
#> 4  1.5952808 0  0.6672881    2 0.6894695
#> 5  0.3295078 1  0.7092303    2 0.5410956
#> 6 -0.8204684 0 -2.0350035    2 0.3988560

# Heterogeneous effect with an explicit covariate spec
data2 <- causalsim(
  n = 500,
  covariates = list(
    W = covar("normal", role = "confounder"),
    V = covar("binary", role = "effect_modifier", prob = 0.4)
  ),
  effect = function(V) 2 + 1.5 * V,
  propensity = function(W) plogis(0.5 * W),
  seed = 42L
)
head(data2)
#>            W V A         Y .tau        .p
#> 1  1.3709584 1 1 5.8250585  3.5 0.6649605
#> 2 -0.5646982 0 1 2.5241222  2.0 0.4298780
#> 3  0.3631284 1 1 4.4707334  3.5 0.5452668
#> 4  0.6328626 0 0 0.3769734  2.0 0.5784543
#> 5  0.4042683 0 1 1.0040666  2.0 0.5503622
#> 6 -0.1061245 0 1 1.4025171  2.0 0.4867375
```
