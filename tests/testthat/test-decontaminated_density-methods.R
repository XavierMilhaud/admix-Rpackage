test_that("print.decontaminated_density displays information", {
  obj <- decontaminated_density(sample1 = rnorm(100), admixMod = admix_model("norm", list(mean = 0, sd = 1)), estim.p = 0.8)
  expect_output(print(obj), "Statistics about the estimated decontaminated density function")
})


test_that("summary.decontaminated_density displays support information", {
  obj <- decontaminated_density(sample1 = rnorm(100), admixMod = admix_model("norm", list(mean = 0, sd = 1)), estim.p = 0.8)
  expect_output(summary(obj), "Type of support: Continuous")
})


test_that("plot.decontaminated_density works for continuous data", {
  obj <- decontaminated_density(sample1 = rnorm(100), admixMod = admix_model("norm", list(mean = 0, sd = 1)), estim.p = 0.8)
  expect_no_error(plot(obj))
})


test_that("plot.decontaminated_density works for discrete data", {
  obj <- decontaminated_density(sample1 = rpois(100, 2), admixMod = admix_model("pois",list(lambda = 1)), estim.p = 0.8)
  expect_no_error(plot(obj))
})
