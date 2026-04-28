# causalsim 0.1.0

Initial CRAN release.

- `causalsim_dgp()` for defining causal data generating processes with known ground truth
- `causalsim_draw()` for simulating individual datasets from a DGP
- `causalsim_eval()` for evaluating estimator performance over repeated replications
- `causalsim_grid()` for sweeping evaluation across a Cartesian product of DGP parameters
- `covar()` for specifying covariate distributions and causal roles
- `summary()` and `plot()` methods for `causalsim_eval` results
