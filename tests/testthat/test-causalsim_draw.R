# ── Return type and shape ─────────────────────────────────────────────────────

test_that("causalsim_draw returns a data frame", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  expect_s3_class(causalsim_draw(dgp), "data.frame")
})

test_that("causalsim_draw returns dgp$n rows", {
  dgp <- causalsim_dgp(n = 250, n_confounders = 1, effect = 1)
  expect_equal(nrow(causalsim_draw(dgp)), 250L)
})

test_that("output has covariate, treatment, outcome, and meta columns", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  d   <- causalsim_draw(dgp)
  expect_true(all(c("W", "A", "Y", ".tau", ".p") %in% names(d)))
})

test_that("multi-covariate DGP produces correctly named columns", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 2,
                       n_noise = 1, effect = 1)
  d   <- causalsim_draw(dgp)
  expect_true(all(c("W1", "W2", "X", "A", "Y") %in% names(d)))
})

# ── Column types and ranges ───────────────────────────────────────────────────

test_that("treatment A is binary (0 or 1)", {
  dgp <- causalsim_dgp(n = 200, n_confounders = 1, effect = 1)
  d   <- causalsim_draw(dgp)
  expect_true(all(d$A %in% c(0L, 1L)))
})

test_that(".p values are in [0, 1]", {
  dgp <- causalsim_dgp(n = 200, n_confounders = 1, effect = 1)
  d   <- causalsim_draw(dgp)
  expect_true(all(d$.p >= 0 & d$.p <= 1))
})

test_that("Y is numeric", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  expect_type(causalsim_draw(dgp)$Y, "double")
})

# ── Reproducibility ───────────────────────────────────────────────────────────

test_that("same seed produces identical draws", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  d1  <- causalsim_draw(dgp, seed = 42)
  d2  <- causalsim_draw(dgp, seed = 42)
  expect_identical(d1, d2)
})

test_that("different seeds produce different draws", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  d1  <- causalsim_draw(dgp, seed = 1)
  d2  <- causalsim_draw(dgp, seed = 2)
  expect_false(isTRUE(all.equal(d1$Y, d2$Y)))
})

# ── Ground truth columns ──────────────────────────────────────────────────────

test_that(".tau is constant for scalar effect", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 3)
  d   <- causalsim_draw(dgp, seed = 1)
  expect_true(all(d$.tau == 3))
})

test_that(".tau varies for heterogeneous effect", {
  dgp <- causalsim_dgp(n = 200, n_confounders = 1,
                       effect = function(W) 2 + W)
  d   <- causalsim_draw(dgp, seed = 1)
  expect_gt(var(d$.tau), 0)
})

test_that("mean .tau is close to true_ate for large draw", {
  # E[2 + W] = 2 when W ~ N(0,1)
  set.seed(7)
  dgp <- causalsim_dgp(n = 10000, n_confounders = 1,
                       effect = function(W) 2 + W)
  d   <- causalsim_draw(dgp)
  expect_equal(mean(d$.tau), dgp$true_ate, tolerance = 0.05)
})

# ── Confounding behaviour ─────────────────────────────────────────────────────

test_that("a default confounder actually induces naive bias", {
  # With default propensity + baseline presets, a confounder must enter both
  # the treatment and outcome models, so a naive difference in means is biased
  # away from the true ATE.
  dgp <- causalsim_dgp(n = 50000, n_confounders = 1, effect = 2)
  d   <- causalsim_draw(dgp, seed = 1)
  naive <- mean(d$Y[d$A == 1]) - mean(d$Y[d$A == 0])
  expect_gt(abs(naive - dgp$true_ate), 0.05)
})

test_that("an RCT (scalar propensity) leaves the naive estimate unbiased", {
  # Randomized treatment: even with a preset baseline the confounder cannot
  # induce bias, because assignment does not depend on it.
  dgp <- causalsim_dgp(n = 50000, n_confounders = 1, effect = 2,
                       propensity = 0.5)
  d   <- causalsim_draw(dgp, seed = 1)
  naive <- mean(d$Y[d$A == 1]) - mean(d$Y[d$A == 0])
  expect_lt(abs(naive - dgp$true_ate), 0.05)
})

# ── No-covariate DGP ─────────────────────────────────────────────────────────

test_that("draw works with no covariates and scalar functions", {
  dgp <- causalsim_dgp(n = 100, effect = 2, propensity = 0.5, baseline = 0)
  d   <- causalsim_draw(dgp, seed = 1)
  expect_equal(nrow(d), 100L)
  expect_true(all(d$.tau == 2))
  expect_true(all(d$.p  == 0.5))
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("non-causalsim_dgp input errors", {
  expect_error(causalsim_draw(list()), "causalsim_dgp")
})
