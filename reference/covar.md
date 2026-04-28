# Define a covariate for a causal DGP

Constructs a fully-specified covariate descriptor for use in
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md).
Every covariate has a distribution family, optional distribution
parameters, and one or more causal roles that control how it enters the
automatically-generated treatment and outcome models.

## Usage

``` r
covar(dist = "normal", role = "confounder", ...)
```

## Arguments

- dist:

  Character. Distribution family. One of `"normal"`, `"binary"`,
  `"uniform"`.

- role:

  Character vector. Causal role(s). One or more of:

  `"confounder"`

  :   Enters both the propensity model and the outcome baseline; creates
      confounding bias in naive estimators.

  `"instrument"`

  :   Enters the propensity model only (exclusion restriction holds).

  `"effect_modifier"`

  :   Available for use in the `effect` function; excluded from the
      propensity model unless also a `"confounder"`.

  `"noise"`

  :   Independent of both treatment and outcome; adds variance without
      confounding.

  Roles are not mutually exclusive:
  `role = c("confounder", "effect_modifier")` specifies a variable that
  both confounds and moderates the effect.

- ...:

  Distribution parameters.

  - `"normal"`: `mean` (default 0), `sd` (default 1)

  - `"binary"`: `prob` (default 0.5)

  - `"uniform"`: `min` (default 0), `max` (default 1)

## Value

An S3 object of class `causalsim_covar`.

## Examples

``` r
# Standard normal confounder
covar("normal", role = "confounder", mean = 0, sd = 1)
#> <causalsim_covar>  dist: normal(mean=0, sd=1)  role: [confounder]

# Binary effect modifier
covar("binary", role = "effect_modifier", prob = 0.4)
#> <causalsim_covar>  dist: binary(prob=0.4)      role: [effect_modifier]

# Variable that both confounds and moderates the effect
covar("normal", role = c("confounder", "effect_modifier"))
#> <causalsim_covar>  dist: normal                role: [confounder, effect_modifier]
```
