test_that("estim_BVdk fails when admixMod is invalid", {
  x <- rnorm(20)
  expect_error(estim_BVdk(samples = x, admixMod = list()),
               "Argument 'admixMod' is not correctly specified")
})


test_that("estim_BVdk uses upper initialization branch", {
  samples <- rnorm(100, mean = 10)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_BVdk(samples = samples, admixMod = admixMod)

  expect_s3_class(res, "estim_BVdk")
})


test_that("estim_BVdk works with Nelder-Mead optimization", {
  samples <- rnorm(50)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_BVdk(samples = samples, admixMod = admixMod, method = "Nelder-Mead")

  expect_equal(res$optim_method, "Nelder-Mead")
})


test_that("print.estim_BVdk prints expected output", {
  skip_on_cran()
  samples <- rnorm(50)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_BVdk(samples = samples, admixMod = admixMod, compute_var = TRUE)

  expect_output(print.estim_BVdk(res), "Estimation \\(BVdk method\\)")
  expect_output(print.estim_BVdk(res), "Mixing weight")
  expect_output(print.estim_BVdk(res), "Optimization method")
  expect_output(print.estim_BVdk(res), "Variance")
})


test_that("summary.estim_BVdk prints summary", {
  samples <- rnorm(50)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_BVdk(samples, admixMod)

  expect_output(summary.estim_BVdk(res), "Known component")
  expect_output(summary.estim_BVdk(res), "Optimization")
})


test_that("summary.estim_BVdk hides call when requested", {
  samples <- rnorm(50)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- estim_BVdk(samples, admixMod)
  txt <- capture.output(summary.estim_BVdk(res, show.call = FALSE))

  expect_false(any(grepl("Call:", txt)))
})


test_that("Donsker_correl_old returns valid value", {
  x <- c(1,2,3,4,5)
  val <- Donsker_correl_old(u = 2, v = 4, obs.data = x)

  expect_true(is.numeric(val))
  expect_length(val, 1)
  expect_gte(val, 0)
})


test_that("BVdk_ML_varCov_estimators returns positive variance", {
  mixt1 <- twoComp_mixt(n = 400, weight = 0.9, comp.dist = list("norm", "norm"),
                         comp.param = list(list("mean" = 4, "sd" = 1), list("mean" = 7, "sd" = 0.5)))
  data1 <- get_mixture_data(mixt1)
  admixMod1 <- admix_model(knownComp_dist = mixt1$comp.dist[[2]], knownComp_param = mixt1$comp.param[[2]])
  res <- BVdk_ML_varCov_estimators(data = data1, admixMod = admixMod1, hat_w = 0.87,
                                   hat_loc = 3.8, hat_var = 1.25)

  expect_true(is.numeric(res))
  expect_length(res, 1)
  expect_gt(res, 0)
})
