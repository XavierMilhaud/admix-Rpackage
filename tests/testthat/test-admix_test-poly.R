test_that("admix_test dispatches to gaussianity_test for one sample", {
  fake_test <- list(statistic = 1.2, p.value = 0.15)
  local_mocked_bindings(
    gaussianity_test = function(sample, admixMod, conf_level, ...) { fake_test }
  )

  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_test(samples = list(rnorm(100)), admixMod = list(mod), test_method = "poly")

  expect_s3_class(res, "gaussianity_test")
  expect_s3_class(res, "htest")
  expect_equal(res$statistic, 1.2)
})


test_that("admix_test dispatches to orthobasis_test for two samples", {
  fake_test <- list(statistic = 2.5, p.value = 0.03)

  local_mocked_bindings(
    orthobasis_test = function(samples, admixMod, conf_level, ...) { fake_test }
  )

  mod1 <- admix_model("norm", list(mean = 0, sd = 1))
  mod2 <- admix_model("norm", list(mean = 1, sd = 1))
  res <- admix_test(samples = list(rnorm(100), rnorm(120)), admixMod = list(mod1, mod2), test_method = "poly")

  expect_s3_class(res, "orthobasis_test")
  expect_s3_class(res, "htest")
  expect_equal(res$p.value, 0.03)
})


test_that("poly testing rejects more than two samples", {
  mod <- admix_model("norm", list(mean = 0, sd = 1))
  expect_error(
    admix_test(samples = list(rnorm(10), rnorm(10), rnorm(10)), admixMod = list(mod, mod, mod), test_method = "poly"),
    "involves at most TWO samples"
  )
})


test_that("poly testing prints information message", {
  fake_test <- list(statistic = 1, p.value = 0.5)
  local_mocked_bindings(
    gaussianity_test = function(sample, admixMod, conf_level, ...) { fake_test }
  )

  mod <- admix_model("norm", list(mean = 0, sd = 1))

  expect_message(
    admix_test(samples = list(rnorm(100)), admixMod = list(mod), test_method = "poly"),
    "Default estimation method is 'BVdk'"
  )
})
