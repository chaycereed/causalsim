# Plot the distribution of estimates from a causalsim_eval result

Draws a histogram of per-replication estimates with vertical lines
marking the true ATE (solid) and the mean estimate (dashed).

## Usage

``` r
# S3 method for class 'causalsim_eval'
plot(x, ...)
```

## Arguments

- x:

  A `causalsim_eval` object.

- ...:

  Additional arguments passed to
  [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html).

## Value

`x`, invisibly.
