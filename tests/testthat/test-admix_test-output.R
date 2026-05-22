test_that("admix_test stores function call", {
  fake_test <- list(statistic = 1, p.value = 0.5)
  local_mocked_bindings(
    gaussianity_test = function(sample, admixMod, conf_level, ...) { fake_test }
  )

  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_test(samples = list(rnorm(10)), admixMod = list(mod), test_method = "poly")

  expect_true(is.call(res$call))
})


test_that("admix_test runs on simple gaussian example", {
  set.seed(123)
  x <- c(rnorm(200, -2, 1), rnorm(100, 0, 1))
  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_test(samples = list(x), admixMod = list(mod), test_method = "poly", conf_level = 0.95)

  expect_s3_class(res, "gaussianity_test")
  expect_true(is.numeric(res$p.value))
})
