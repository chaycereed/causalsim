# Create a causal data generating process

Defines a causal DGP with known ground truth. Covariates can be
specified via shorthand count arguments (Option B), an explicit named
list of
[`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
objects (Option A), or both combined.

## Usage

``` r
causalsim_dgp(
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
  mc_draws = 10000L
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

## Value

An S3 object of class `causalsim_dgp` with components:

- `n`:

  Sample size (integer)

- `covar_spec`:

  Named list of
  [`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
  objects

- `effect_fn`:

  Normalized effect function

- `propensity_fn`:

  Normalized propensity function

- `baseline_fn`:

  Normalized baseline function

- `sigma`:

  Outcome noise standard deviation

- `true_ate`:

  True ATE: exact for scalar, MC approximation otherwise

- `heterogeneous`:

  Logical; `TRUE` if `effect` was a function

- `mc_draws`:

  Monte Carlo draws used (integer)

## Details

### Structural model

    W ~ covariate_spec
    A ~ Bernoulli(propensity(W))
    Y  = baseline(W) + effect(W) * A + N(0, sigma^2)

`A` is used for treatment throughout to avoid collision with R's
built-in `T` alias.

### Function calling convention

The `effect`, `propensity`, and `baseline` functions are called with
named arguments matching covariate names, not a data frame. Write:

    effect = function(W) 2 + 1.5 * W
    propensity = function(W1, W2) plogis(0.3 * W1 + 0.5 * W2)

Every argument name is validated against the DGP's covariate spec at
construction time, so mismatches surface immediately rather than at draw
time. The propensity function is additionally evaluated on a small test
draw to confirm it returns values in \[0, 1\].

### True ATE

For scalar `effect` the true ATE is exact. For function `effect` it is
approximated via Monte Carlo over `mc_draws` draws from the covariate
distribution.

## Examples

``` r
# Minimal: one confounder, constant effect, moderate confounding
dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2)
dgp
#> <causalsim_dgp>
#>   n            : 500
#>   true ATE     : 2.0000
#>   heterogeneous: FALSE
#>   sigma        : 1.00
#>   covariates   :
#>     W       normal  [confounder]

# Heterogeneous effect, explicit covariate spec
dgp2 <- causalsim_dgp(
  n = 500,
  covariates = list(
    W = covar("normal", role = "confounder"),
    V = covar("binary", role = "effect_modifier", prob = 0.4)
  ),
  effect = function(V) 2 + 1.5 * V,
  propensity = function(W) plogis(0.5 * W),
  baseline = function(W) 1.5 * W
)

# Mixed: shorthand confounders + explicit instrument + RCT propensity
dgp3 <- causalsim_dgp(
  n = 1000,
  n_confounders = 2,
  covariates = list(Z = covar("normal", role = "instrument")),
  effect = 1,
  propensity = 0.5
)
```
