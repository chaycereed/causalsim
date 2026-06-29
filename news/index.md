# Changelog

## causalsim 0.1.0

Initial CRAN release.

- [`causalsim()`](https://chaycereed.github.io/causalsim/reference/causalsim.md)
  for simulating a dataset in a single call (convenience wrapper around
  [`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
  and
  [`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md))
- [`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
  for defining causal data generating processes with known ground truth
- [`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
  for simulating individual datasets from a DGP
- [`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
  for evaluating estimator performance over repeated replications
- [`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
  for sweeping evaluation across a Cartesian product of DGP parameters
- [`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
  for specifying covariate distributions and causal roles
- [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  `causalsim_eval` results
