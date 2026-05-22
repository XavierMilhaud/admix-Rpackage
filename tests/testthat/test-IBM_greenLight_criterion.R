test_that("IBM_greenLight_criterion requires exactly 2 samples", {
  estim_obj <- list(integ.supp = 1:10, estimated_mixing_weights = c(0.5, 0.6), p.X.fixed = 0.4)

  expect_error(
    IBM_greenLight_criterion(estim_obj = estim_obj, samples = list(rnorm(10)), admixMod = list(), alpha = 0.05),
    "exactly TWO")
})


test_that("green light criterion returns expected structure", {
  skip_if_not_installed("mockery")
  estim_obj <- list(integ.supp = seq(0, 1, length.out = 20), estimated_mixing_weights = c(0.4, 0.5), p.X.fixed = 0.3)
  samples <- list(rnorm(100), rnorm(120))
  admixMod <- list(structure(list(), class = "admix_model"), structure(list(), class = "admix_model"))
  mockery::stub(
    IBM_greenLight_criterion, "IBM_estimVarCov_gaussVect", function(...) diag(c(0.01, 0.01))
  )
  res <- IBM_greenLight_criterion(estim_obj = estim_obj, samples = samples, admixMod = admixMod)

  expect_type(res, "list")
  expect_named(res, c("green_light", "conf_interval_p1", "conf_interval_p2"))
  expect_type(res$green_light, "logical")
  expect_length(res$conf_interval_p1, 2)
  expect_length(res$conf_interval_p2, 2)
})


test_that("green light criterion can reject impossible weights", {
  skip_if_not_installed("mockery")
  estim_obj <- list(integ.supp = seq(0, 1, length.out = 20), estimated_mixing_weights = c(1.3, 1.4), p.X.fixed = 0.3)
  samples <- list(rnorm(100), rnorm(100))
  admixMod <- list(structure(list(), class = "admix_model"), structure(list(), class = "admix_model"))
  mockery::stub(
    IBM_greenLight_criterion, "IBM_estimVarCov_gaussVect", function(...) diag(c(0.001, 0.001))
  )
  res <- IBM_greenLight_criterion(estim_obj, samples, admixMod)

  expect_false(res$green_light)
})


test_that("green light criterion handles one estimated weight", {
  skip_if_not_installed("mockery")
  estim_obj <- list(integ.supp = seq(0, 1, length.out = 20), estimated_mixing_weights = 0.7, p.X.fixed = 0.4)
  samples <- list(rnorm(100), rnorm(100))
  admixMod <- list(structure(list(), class = "admix_model"), structure(list(), class = "admix_model"))
  mockery::stub(
    IBM_greenLight_criterion, "IBM_estimVarCov_gaussVect", function(...) matrix(0.01, 1, 1))
  res <- IBM_greenLight_criterion(estim_obj, samples, admixMod)

  expect_null(res$conf_interval_p1)
  expect_length(res$conf_interval_p2, 2)
})
