test_that("admix_model creates a valid object for normal distribution", {
  obj <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  expect_s3_class(obj, "admix_model")
  expect_equal(obj$comp.dist$known, "norm")
  expect_equal(obj$comp.param$known, list(mean = 0, sd = 1))
  expect_null(obj$comp.dist$unknown)
  expect_null(obj$comp.param$unknown)
})

test_that("admix_model rejects invalid knownComp_dist", {
  expect_error(admix_model(knownComp_dist = 1, knownComp_param = list(mean = 0, sd = 1)),
               "`knownComp_dist` must be a character string of length 1.")
  expect_error(admix_model(knownComp_dist = c("norm", "exp"), knownComp_param = list(mean = 0, sd = 1)),
               "`knownComp_dist` must be a character string of length 1.")
})

test_that("admix_model rejects invalid parameter lists", {
  expect_error(admix_model(knownComp_dist = "norm", knownComp_param = c(mean = 0, sd = 1)),
               "`knownComp_param` must be a named list.")
  expect_error(admix_model(knownComp_dist = "norm", knownComp_param = list(0, 1)),
               "`knownComp_param` must be a named list.")
})

test_that("admix_model works with several distributions", {
  obj_exp <- admix_model("exp", list(rate = 2))
  expect_equal(obj_exp$comp.dist$known, "exp")
  obj_pois <- admix_model("pois", list(lambda = 5))
  expect_equal(obj_pois$comp.dist$known, "pois")
})

