test_that("covar returns correct S3 class", {
  cv <- covar("normal", role = "confounder")
  expect_s3_class(cv, "causalsim_covar")
})

test_that("covar stores dist, role, and params", {
  cv <- covar("binary", role = "effect_modifier", prob = 0.3)
  expect_equal(cv$dist, "binary")
  expect_equal(cv$role, "effect_modifier")
  expect_equal(cv$params$prob, 0.3)
})

test_that("covar accepts multiple roles", {
  cv <- covar("normal", role = c("confounder", "effect_modifier"))
  expect_equal(cv$role, c("confounder", "effect_modifier"))
})

test_that("covar rejects invalid dist", {
  expect_error(covar("gamma"), "should be one of")
})

test_that("covar rejects invalid role", {
  expect_error(covar("normal", role = "collider"), "should be one of")
})

test_that("print.causalsim_covar runs without error", {
  cv <- covar("normal", role = "confounder", mean = 0, sd = 1)
  expect_output(print(cv), "causalsim_covar")
  expect_invisible(print(cv))
})
