# Package index

## Data Generating Processes

Define and draw from causal DGPs with known ground truth.

- [`causalsim_dgp()`](https://chaycereed.github.io/causalsim/reference/causalsim_dgp.md)
  : Create a causal data generating process
- [`covar()`](https://chaycereed.github.io/causalsim/reference/covar.md)
  : Define a covariate for a causal DGP
- [`causalsim_draw()`](https://chaycereed.github.io/causalsim/reference/causalsim_draw.md)
  : Draw a dataset from a causal DGP

## Evaluation

Measure estimator performance against the known ground truth.

- [`causalsim_eval()`](https://chaycereed.github.io/causalsim/reference/causalsim_eval.md)
  : Evaluate a causal estimator against a known DGP
- [`summary(`*`<causalsim_eval>`*`)`](https://chaycereed.github.io/causalsim/reference/summary.causalsim_eval.md)
  : Summarise a causalsim_eval result
- [`plot(`*`<causalsim_eval>`*`)`](https://chaycereed.github.io/causalsim/reference/plot.causalsim_eval.md)
  : Plot the distribution of estimates from a causalsim_eval result
- [`causalsim_grid()`](https://chaycereed.github.io/causalsim/reference/causalsim_grid.md)
  : Evaluate an estimator across a parameter grid
