# Summarise a causalsim_eval result

Returns a structured summary of estimator performance including the
metrics table and the distribution of per-replication estimates (mean,
SD, median, and 10th / 90th percentiles).

## Usage

``` r
# S3 method for class 'causalsim_eval'
summary(object, ...)
```

## Arguments

- object:

  A `causalsim_eval` object.

- ...:

  Ignored.

## Value

An S3 object of class `causalsim_eval_summary` with components:

- `metrics`:

  The metrics data frame from the original result.

- `mean_estimate`:

  Mean of per-replication estimates.

- `sd_estimate`:

  Standard deviation of per-replication estimates.

- `median_estimate`:

  Median of per-replication estimates.

- `p10_estimate`:

  10th percentile of per-replication estimates.

- `p90_estimate`:

  90th percentile of per-replication estimates.

- `true_ate`:

  True ATE from the DGP.

- `reps`:

  Number of replications.
