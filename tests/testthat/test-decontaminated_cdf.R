test_that("decontaminated_cdf validates admix model", {
  expect_error(
    decontaminated_cdf(sample1 = rnorm(10), admixMod = "bad", estim.p = 0.8),
    "not correctly specified")
})


test_that("decontaminated_cdf returns a stepfun for continuous data", {
  obj <- decontaminated_cdf(sample1 = rnorm(200), admixMod = admix_model("norm",list(mean = 0, sd = 1)), estim.p = 0.8)
  expect_s3_class(obj, "stepfun")

  vals <- obj(c(-1,0,1))
  expect_length(vals, 3)
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1))
})


test_that("decontaminated_cdf works for discrete data", {
  obj <- decontaminated_cdf(sample1 = rpois(200, 2), admixMod = admix_model("pois",list(lambda = 1)), estim.p = 0.7)
  vals <- obj(c(0,1,2))

  expect_length(vals, 3)
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1))
})
