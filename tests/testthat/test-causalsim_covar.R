test_that("causalsim_covar returns correct S3 class", {
  cv <- causalsim_covar("normal", role = "confounder")
  expect_s3_class(cv, "causalsim_covar")
})

test_that("causalsim_covar stores dist, role, and params", {
  cv <- causalsim_covar("binary", role = "effect_modifier", prob = 0.3)
  expect_equal(cv$dist, "binary")
  expect_equal(cv$role, "effect_modifier")
  expect_equal(cv$params$prob, 0.3)
})

test_that("causalsim_covar accepts multiple roles", {
  cv <- causalsim_covar("normal", role = c("confounder", "effect_modifier"))
  expect_equal(cv$role, c("confounder", "effect_modifier"))
})

test_that("causalsim_covar rejects invalid dist", {
  expect_error(causalsim_covar("gamma"), "should be one of")
})

test_that("causalsim_covar rejects invalid role", {
  expect_error(causalsim_covar("normal", role = "collider"), "should be one of")
})

test_that("print.causalsim_covar runs without error", {
  cv <- causalsim_covar("normal", role = "confounder", mean = 0, sd = 1)
  expect_output(print(cv), "causalsim_covar")
  expect_invisible(print(cv))
})
