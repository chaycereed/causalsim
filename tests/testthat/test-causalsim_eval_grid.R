# ── Shared fixtures ───────────────────────────────────────────────────────────

dgp_base <- causalsim_dgp(
  n             = 100,
  n_confounders = 1,
  effect        = 2,
  propensity    = 0.5
)

# Constant estimator: deterministic, CI always brackets true ATE = 2
const_est <- function(data) {
  c(estimate = 2, ci_lower = 1, ci_upper = 3)
}

# ── Return type and structure ─────────────────────────────────────────────────

test_that("causalsim_grid returns correct S3 class", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L)),
                      reps = 5L, metrics = "bias", seed = 1L)
  expect_s3_class(g, "causalsim_grid")
})

test_that("results has vary columns + metric + value + se", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L)),
                      reps = 5L, metrics = "bias", seed = 1L)
  expect_true(all(c("n", "metric", "value", "se") %in% names(g$results)))
})

test_that("results has n_cells x n_metrics rows", {
  g <- causalsim_grid(
    dgp_base, const_est,
    vary    = list(n = c(100L, 200L), sigma = c(0.5, 1.0)),
    reps    = 5L,
    metrics = c("bias", "rmse"),
    seed    = 1L
  )
  expect_equal(nrow(g$results), 8L)  # 2 x 2 cells x 2 metrics
})

test_that("grid parameter values appear correctly in results", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L, 300L)),
                      reps = 5L, metrics = "bias", seed = 1L)
  expect_setequal(unique(g$results$n), c(100L, 200L, 300L))
})

test_that("vary parameter names appear in $vary", {
  g <- causalsim_grid(
    dgp_base, const_est,
    vary = list(n = c(100L, 200L), sigma = c(0.5, 1.0)),
    reps = 5L, metrics = "bias", seed = 1L
  )
  expect_equal(g$vary, c("n", "sigma"))
})

test_that("$grid has one row per cell", {
  g <- causalsim_grid(
    dgp_base, const_est,
    vary = list(n = c(100L, 200L), sigma = c(0.5, 1.0)),
    reps = 5L, metrics = "bias", seed = 1L
  )
  expect_equal(nrow(g$grid), 4L)
})

test_that("reps is stored on the result", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L)),
                      reps = 7L, metrics = "bias", seed = 1L)
  expect_equal(g$reps, 7L)
})

# ── Metric values ─────────────────────────────────────────────────────────────

test_that("bias is 0 for constant-exact estimator across all cells", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L)),
                      reps = 5L, metrics = "bias", seed = 1L)
  bias_vals <- g$results$value[g$results$metric == "bias"]
  expect_true(all(bias_vals == 0))
})

# ── Reproducibility ───────────────────────────────────────────────────────────

test_that("same seed gives identical results", {
  g1 <- causalsim_grid(dgp_base, const_est,
                       vary = list(n = c(100L, 200L)),
                       reps = 5L, metrics = "bias", seed = 42L)
  g2 <- causalsim_grid(dgp_base, const_est,
                       vary = list(n = c(100L, 200L)),
                       reps = 5L, metrics = "bias", seed = 42L)
  expect_identical(g1$results, g2$results)
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("non-causalsim_dgp input errors", {
  expect_error(
    causalsim_grid(list(), const_est, vary = list(n = 100L), reps = 5L),
    "causalsim_dgp"
  )
})

test_that("non-function estimator errors", {
  expect_error(
    causalsim_grid(dgp_base, "ols", vary = list(n = 100L), reps = 5L),
    "`estimator`"
  )
})

test_that("empty vary errors", {
  expect_error(
    causalsim_grid(dgp_base, const_est, vary = list(), reps = 5L),
    "non-empty"
  )
})

test_that("unnamed vary errors", {
  expect_error(
    causalsim_grid(dgp_base, const_est,
                   vary = list(c(100L, 200L)), reps = 5L),
    "named"
  )
})

test_that("invalid vary parameter name errors", {
  expect_error(
    causalsim_grid(dgp_base, const_est,
                   vary = list(foo = c(1, 2)), reps = 5L),
    "not a valid"
  )
})

test_that("non-atomic vary element errors", {
  expect_error(
    causalsim_grid(
      dgp_base, const_est,
      vary = list(effect = list(function(W) W, function(W) 2 * W)),
      reps = 5L
    ),
    "atomic"
  )
})

# ── Print method ──────────────────────────────────────────────────────────────

test_that("print runs without error and returns invisibly", {
  g <- causalsim_grid(dgp_base, const_est,
                      vary = list(n = c(100L, 200L)),
                      reps = 5L, metrics = "bias", seed = 1L)
  expect_invisible(print(g))
  expect_output(print(g), "causalsim_grid")
  expect_output(print(g), "vary")
})
