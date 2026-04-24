# ── Shared fixtures ───────────────────────────────────────────────────────────

# DGP: RCT (propensity = 0.5) so OLS is unbiased; true ATE = 2
dgp_rct <- causalsim_dgp(
  n             = 300,
  n_confounders = 1,
  effect        = 2,
  propensity    = 0.5
)

# Constant estimator: always returns exactly the true ATE with a wide CI
const_estimator <- function(data) {
  c(estimate = 2, ci_lower = 1, ci_upper = 3)
}

# Biased estimator: always overshoots by 0.5
biased_estimator <- function(data) {
  c(estimate = 2.5, ci_lower = 2, ci_upper = 3)
}

# Estimator returning a data frame instead of a vector
df_estimator <- function(data) {
  data.frame(estimate = 2, ci_lower = 1, ci_upper = 3)
}

# Estimator with no CI (bias/rmse only)
noci_estimator <- function(data) {
  c(estimate = 2)
}

# ── Return type and structure ─────────────────────────────────────────────────

test_that("causalsim_eval returns causalsim_eval class", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L,
                           seed = 1L)
  expect_s3_class(result, "causalsim_eval")
})

test_that("result$metrics has metric, value, se columns", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L, seed = 1L)
  expect_named(result$metrics, c("metric", "value", "se"))
})

test_that("result$draws has one row per rep", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 1L)
  expect_equal(nrow(result$draws), 20L)
})

test_that("result$draws has estimate column", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L, seed = 1L)
  expect_true("estimate" %in% names(result$draws))
})

test_that("result$reps matches requested reps", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 15L, seed = 1L)
  expect_equal(result$reps, 15L)
})

test_that("result$true_ate matches dgp$true_ate", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L, seed = 1L)
  expect_equal(result$true_ate, dgp_rct$true_ate)
})

# ── Metric values ─────────────────────────────────────────────────────────────

test_that("bias is 0 for constant-exact estimator", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 50L,
                           metrics = "bias", seed = 1L)
  expect_equal(result$metrics$value[result$metrics$metric == "bias"], 0)
})

test_that("bias is approximately 0.5 for consistently biased estimator", {
  result <- causalsim_eval(dgp_rct, biased_estimator, reps = 50L,
                           metrics = "bias", seed = 1L)
  expect_equal(
    result$metrics$value[result$metrics$metric == "bias"], 0.5,
    tolerance = 1e-10
  )
})

test_that("rmse is 0 for constant-exact estimator", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 50L,
                           metrics = "rmse", seed = 1L)
  expect_equal(result$metrics$value[result$metrics$metric == "rmse"], 0)
})

test_that("coverage is 1 when CI always brackets true_ate", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 50L,
                           metrics = "coverage", seed = 1L)
  expect_equal(
    result$metrics$value[result$metrics$metric == "coverage"], 1
  )
})

test_that("power is 1 when CI always excludes zero", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 50L,
                           metrics = "power", seed = 1L)
  expect_equal(
    result$metrics$value[result$metrics$metric == "power"], 1
  )
})

test_that("only requested metrics are returned", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L,
                           metrics = c("bias", "rmse"), seed = 1L)
  expect_equal(sort(result$metrics$metric), c("bias", "rmse"))
})

# ── Data frame estimator ──────────────────────────────────────────────────────

test_that("estimator returning a data frame is accepted", {
  expect_no_error(
    causalsim_eval(dgp_rct, df_estimator, reps = 10L,
                   metrics = "bias", seed = 1L)
  )
})

# ── Named list estimator ──────────────────────────────────────────────────────

test_that("estimator returning a named numeric list is accepted", {
  list_estimator <- function(data) {
    list(estimate = 2, ci_lower = 1, ci_upper = 3)
  }
  expect_no_error(
    causalsim_eval(dgp_rct, list_estimator, reps = 10L,
                   metrics = "bias", seed = 1L)
  )
})

test_that("list estimator gives same result as equivalent vector estimator", {
  vec_est  <- function(data) c(estimate = 2, ci_lower = 1, ci_upper = 3)
  list_est <- function(data) list(estimate = 2, ci_lower = 1, ci_upper = 3)
  r1 <- causalsim_eval(dgp_rct, vec_est,  reps = 10L,
                       metrics = "bias", seed = 5L)
  r2 <- causalsim_eval(dgp_rct, list_est, reps = 10L,
                       metrics = "bias", seed = 5L)
  expect_equal(r1$metrics, r2$metrics)
})

test_that("inner names on list values do not override outer names", {
  inner_named_est <- function(data) {
    fit <- lm(Y ~ A, data = data)
    list(
      estimate = coef(fit)["A"],
      ci_lower = confint(fit)["A", 1],
      ci_upper = confint(fit)["A", 2]
    )
  }
  result <- causalsim_eval(dgp_rct, inner_named_est, reps = 5L,
                           metrics = "bias", seed = 1L)
  expect_true("bias" %in% result$metrics$metric)
})

test_that("list estimator with non-numeric element errors", {
  bad_list_est <- function(data) {
    list(estimate = 2, ci_lower = "1", ci_upper = 3)
  }
  expect_error(
    causalsim_eval(dgp_rct, bad_list_est, reps = 5L, seed = 1L),
    "not all elements are named numeric"
  )
})

# ── No-CI estimator ───────────────────────────────────────────────────────────

test_that("bias and rmse work with no-CI estimator", {
  expect_no_error(
    causalsim_eval(dgp_rct, noci_estimator, reps = 10L,
                   metrics = c("bias", "rmse"), seed = 1L)
  )
})

test_that("requesting coverage with no-CI estimator errors", {
  expect_error(
    causalsim_eval(dgp_rct, noci_estimator, reps = 10L,
                   metrics = "coverage", seed = 1L),
    "ci_lower"
  )
})

# ── Reproducibility ───────────────────────────────────────────────────────────

test_that("same seed gives identical results", {
  r1 <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 7L)
  r2 <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 7L)
  expect_identical(r1$draws, r2$draws)
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("non-causalsim_dgp input errors", {
  expect_error(causalsim_eval(list(), const_estimator), "causalsim_dgp")
})

test_that("non-function estimator errors", {
  expect_error(causalsim_eval(dgp_rct, "lm"), "`estimator`")
})

test_that("estimator missing 'estimate' field errors", {
  bad_estimator <- function(data) c(coef = 2, lower = 1, upper = 3)
  expect_error(
    causalsim_eval(dgp_rct, bad_estimator, reps = 5L, seed = 1L),
    "'estimate'"
  )
})

test_that("invalid metric name errors", {
  expect_error(
    causalsim_eval(dgp_rct, const_estimator, reps = 5L,
                   metrics = "variance"),
    "should be one of"
  )
})

# ── Print method ──────────────────────────────────────────────────────────────

test_that("print method runs without error and returns invisibly", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L, seed = 1L)
  expect_invisible(print(result))
  expect_output(print(result), "causalsim_eval")
  expect_output(print(result), "true ATE")
})

# ── Summary method ────────────────────────────────────────────────────────────

test_that("summary returns causalsim_eval_summary class", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 1L)
  expect_s3_class(summary(result), "causalsim_eval_summary")
})

test_that("summary has expected fields", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 1L)
  s <- summary(result)
  expect_true(all(c("mean_estimate", "sd_estimate", "median_estimate",
                    "p10_estimate", "p90_estimate",
                    "true_ate", "reps", "metrics") %in% names(s)))
})

test_that("summary mean_estimate is correct for constant estimator", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 1L)
  s <- summary(result)
  expect_equal(s$mean_estimate, 2)
})

test_that("summary prints without error and returns invisibly", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 10L, seed = 1L)
  s <- summary(result)
  expect_invisible(print(s))
  expect_output(print(s), "mean")
  expect_output(print(s), "true ATE")
})

# ── Plot method ───────────────────────────────────────────────────────────────

test_that("plot runs without error and returns invisibly", {
  result <- causalsim_eval(dgp_rct, const_estimator, reps = 20L, seed = 1L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(result))
})
