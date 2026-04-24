# ── Basic construction ────────────────────────────────────────────────────────

test_that("causalsim_dgp returns correct S3 class", {
  dgp <- causalsim_dgp(n = 100, effect = 1)
  expect_s3_class(dgp, "causalsim_dgp")
})

test_that("scalar effect is normalized to a function", {
  dgp <- causalsim_dgp(n = 100, effect = 3)
  expect_true(is.function(dgp$effect_fn))
})

test_that("heterogeneous is FALSE for scalar effect", {
  dgp <- causalsim_dgp(n = 100, effect = 2)
  expect_false(dgp$heterogeneous)
})

test_that("heterogeneous is TRUE for function effect", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = function(W) 1 + W)
  expect_true(dgp$heterogeneous)
})

test_that("n is stored as integer", {
  dgp <- causalsim_dgp(n = 200.0, effect = 1)
  expect_type(dgp$n, "integer")
  expect_equal(dgp$n, 200L)
})

test_that("mc_draws is stored as integer", {
  dgp <- causalsim_dgp(n = 100, effect = 1, mc_draws = 5000)
  expect_type(dgp$mc_draws, "integer")
})

# ── True ATE ──────────────────────────────────────────────────────────────────

test_that("true_ate is exact for scalar effect", {
  dgp <- causalsim_dgp(n = 100, effect = 2.5)
  expect_equal(dgp$true_ate, 2.5)
})

test_that("true_ate is approx correct for zero-mean heterogeneous effect", {
  # E[2 + 1.5 * W] = 2 when W ~ N(0, 1)
  set.seed(42)
  dgp <- causalsim_dgp(
    n            = 100,
    n_confounders = 1,
    effect       = function(W) 2 + 1.5 * W,
    mc_draws     = 100000L
  )
  expect_equal(dgp$true_ate, 2, tolerance = 0.02)
})

# ── sigma ─────────────────────────────────────────────────────────────────────

test_that("sigma is stored on DGP object", {
  dgp <- causalsim_dgp(n = 100, effect = 1, sigma = 2.5)
  expect_equal(dgp$sigma, 2.5)
})

test_that("negative sigma throws error", {
  expect_error(causalsim_dgp(n = 100, effect = 1, sigma = -1), "`sigma`")
})

test_that("zero sigma throws error", {
  expect_error(causalsim_dgp(n = 100, effect = 1, sigma = 0), "`sigma`")
})

# ── propensity ────────────────────────────────────────────────────────────────

test_that("propensity preset string is normalized to a function", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1,
                       propensity = "moderate", effect = 1)
  expect_true(is.function(dgp$propensity_fn))
})

test_that("propensity scalar is normalized to a function", {
  dgp <- causalsim_dgp(n = 100, propensity = 0.3, effect = 1)
  expect_true(is.function(dgp$propensity_fn))
})

test_that("propensity function returning value > 1 errors at construction", {
  expect_error(
    causalsim_dgp(
      n = 100, n_confounders = 1,
      propensity = function(W) 1.5,
      effect     = 1
    ),
    "\\[0, 1\\]"
  )
})

test_that("propensity function returning negative value errors", {
  expect_error(
    causalsim_dgp(
      n = 100, n_confounders = 1,
      propensity = function(W) -0.1,
      effect     = 1
    ),
    "\\[0, 1\\]"
  )
})

test_that("propensity referencing undefined covariate errors", {
  expect_error(
    causalsim_dgp(
      n = 100, n_confounders = 1,
      propensity = function(W, X) plogis(W + X),
      effect     = 1
    ),
    "not defined in the DGP"
  )
})

# ── baseline ──────────────────────────────────────────────────────────────────

test_that("baseline scalar is normalized to a function", {
  dgp <- causalsim_dgp(n = 100, baseline = 2, effect = 1)
  expect_true(is.function(dgp$baseline_fn))
})

test_that("baseline preset string is normalized to a function", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1,
                       baseline = "high", effect = 1)
  expect_true(is.function(dgp$baseline_fn))
})

# ── Covariate spec: shorthand (Option B) ──────────────────────────────────────

test_that("n_confounders = 1 generates covariate named W", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  expect_named(dgp$covar_spec, "W")
  expect_equal(dgp$covar_spec$W$role, "confounder")
})

test_that("n_confounders = 2 generates W1 and W2", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 2, effect = 1)
  expect_named(dgp$covar_spec, c("W1", "W2"))
})

test_that("n_instruments = 1 generates covariate named Z", {
  dgp <- causalsim_dgp(n = 100, n_instruments = 1, effect = 1)
  expect_named(dgp$covar_spec, "Z")
  expect_equal(dgp$covar_spec$Z$role, "instrument")
})

test_that("n_effect_modifiers = 1 generates covariate named V", {
  dgp <- causalsim_dgp(n = 100, n_effect_modifiers = 1, effect = 1)
  expect_named(dgp$covar_spec, "V")
})

# ── Covariate spec: explicit list (Option A) ──────────────────────────────────

test_that("explicit covariates list is stored on DGP", {
  dgp <- causalsim_dgp(
    n = 100,
    covariates = list(
      W = covar("normal", role = "confounder"),
      V = covar("binary", role = "effect_modifier", prob = 0.4)
    ),
    effect = 1
  )
  expect_named(dgp$covar_spec, c("W", "V"))
  expect_s3_class(dgp$covar_spec$W, "causalsim_covar")
})

# ── Covariate spec: mixed (Option C) ─────────────────────────────────────────

test_that("shorthand and explicit covariates merge correctly", {
  dgp <- causalsim_dgp(
    n = 100,
    n_confounders = 1,
    covariates    = list(Z = covar("normal", role = "instrument")),
    effect        = 1
  )
  expect_named(dgp$covar_spec, c("W", "Z"))
})

test_that("name collision between auto and explicit covariates errors", {
  expect_error(
    causalsim_dgp(
      n = 100,
      n_confounders = 1,
      covariates    = list(W = covar("binary", role = "confounder")),
      effect        = 1
    ),
    "appear in both"
  )
})

# ── Effect function validation ────────────────────────────────────────────────

test_that("effect referencing undefined covariate errors at construction", {
  expect_error(
    causalsim_dgp(n = 100, n_confounders = 1,
                  effect = function(W, X) W + X),
    "not defined in the DGP"
  )
})

test_that("effect referencing defined covariates does not error", {
  expect_no_error(
    causalsim_dgp(n = 100, n_confounders = 1, effect = function(W) 2 + W)
  )
})

# ── Input validation ──────────────────────────────────────────────────────────

test_that("negative n throws error", {
  expect_error(causalsim_dgp(n = -1, effect = 1), "`n`")
})

test_that("zero n throws error", {
  expect_error(causalsim_dgp(n = 0, effect = 1), "`n`")
})

test_that("non-numeric n throws error", {
  expect_error(causalsim_dgp(n = "100", effect = 1), "`n`")
})

test_that("vector effect throws error", {
  expect_error(causalsim_dgp(n = 100, effect = c(1, 2)), "`effect`")
})

test_that("non-numeric non-function effect throws error", {
  expect_error(causalsim_dgp(n = 100, effect = "large"), "preset")
})

test_that("unnamed covariates list throws error", {
  expect_error(
    causalsim_dgp(n = 100, covariates = list(covar("normal"))),
    "must be named"
  )
})

test_that("non-covar element in covariates list throws error", {
  expect_error(
    causalsim_dgp(n = 100, covariates = list(W = list(dist = "normal"))),
    "covar\\(\\) objects"
  )
})

# ── Print method ──────────────────────────────────────────────────────────────

test_that("print method runs without error and returns invisibly", {
  dgp <- causalsim_dgp(n = 100, n_confounders = 1, effect = 1)
  expect_invisible(print(dgp))
  expect_output(print(dgp), "causalsim_dgp")
  expect_output(print(dgp), "true ATE")
  expect_output(print(dgp), "W")
  expect_output(print(dgp), "sigma")
})
