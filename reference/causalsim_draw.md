# Draw a dataset from a causal DGP

Simulates one dataset from a
[`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
object using the structural model:

## Usage

``` r
causalsim_draw(dgp, seed = NULL)
```

## Arguments

- dgp:

  A `causalsim_dgp` object created by
  [`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md).

- seed:

  Integer or `NULL`. If non-null, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before generating
  any random values, making the draw reproducible. Default `NULL`.

## Value

A data frame with `dgp$n` rows and the following columns:

- Covariate columns:

  One column per covariate in `dgp$covar_spec`, named to match the spec
  (e.g. `W`, `W1`, `Z`).

- `A`:

  Binary treatment indicator (0/1).

- `Y`:

  Observed outcome.

- `.tau`:

  Individual treatment effect (CATE). Ground truth.

- `.p`:

  Individual propensity score. Ground truth.

The `.tau` and `.p` columns carry individual-level ground truth and are
prefixed with `.` to distinguish them from observed variables.
[`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
uses them directly without re-estimation.

## Details

    W ~ covariate_spec
    A ~ Bernoulli(propensity(W))
    Y  = baseline(W) + effect(W) * A + N(0, sigma^2)

## Examples

``` r
dgp <- causalsim_dgp(n = 500, n_confounders = 1, effect = 2)
d   <- causalsim_draw(dgp, seed = 1)
head(d)
#>            W A          Y .tau        .p
#> 1 -0.6264538 0 -1.4546914    2 0.4223273
#> 2  0.1836433 0 -0.8456543    2 0.5229393
#> 3 -0.8356286 0 -1.2504797    2 0.3970399
#> 4  1.5952808 0  0.6672881    2 0.6894695
#> 5  0.3295078 1  0.7092303    2 0.5410956
#> 6 -0.8204684 0 -2.0350035    2 0.3988560

# Reproducible draws
d1 <- causalsim_draw(dgp, seed = 42)
d2 <- causalsim_draw(dgp, seed = 42)
identical(d1, d2)  # TRUE
#> [1] TRUE
```
