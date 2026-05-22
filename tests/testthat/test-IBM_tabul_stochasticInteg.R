test_that("IBM_tabul_stochasticInteg requires two samples", {
  expect_error(
    IBM_tabul_stochasticInteg(samples = list(rnorm(10)), admixMod = list()),
    "Must be 2")
})


test_that("IBM_tabul_stochasticInteg returns expected object", {
  skip_on_cran()
  skip_if_not_installed("mockery")

  samples <- list(rnorm(50), rnorm(60))
  admixMod <- list(admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1)),
                   admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1)))
  fake_estim <- list(estimated_mixing_weights = c(0.4, 0.5), integ.supp = seq(-2, 2, length.out = 20), p.X.fixed = 0.3)

  mockery::stub(IBM_tabul_stochasticInteg, "estim_IBM", function(...) fake_estim)
  mockery::stub(IBM_tabul_stochasticInteg, "IBM_empirical_contrast", function(...) 0.5)
  mockery::stub(IBM_tabul_stochasticInteg, "detect_support_type", function(...) "Continuous")
  mockery::stub(IBM_tabul_stochasticInteg, "IBM_normalization_term", function(...) diag(3))
  mockery::stub(IBM_tabul_stochasticInteg, "estimVarCov_empProcess_Rcpp", function(...) diag(10))
  mockery::stub(IBM_tabul_stochasticInteg, "is_equal_knownComp", function(...) FALSE)
  mockery::stub(IBM_tabul_stochasticInteg, "sim_gaussianProcess", function(...) list(traj1 = rnorm(10)))
  res <- IBM_tabul_stochasticInteg(samples = samples, admixMod = admixMod, n_sim_tab = 2, parallel = FALSE)

  expect_type(res, "list")
  expect_named(res, c("U_sim","estimator","contrast_value","integ.points"))
  expect_true(is.numeric(res$U_sim))
})


test_that("tabulated distribution removes error values", {
  U_sim <- c("Error in foo", "list(a=1)", "0.45", "1.2")
  indexes.toRemove <- which(
    (substr(U_sim, 1, 5) == "Error") |
      (substr(U_sim, 1, 5) == "list(")
  )
  cleaned <- as.numeric(U_sim[-indexes.toRemove])

  expect_equal(cleaned, c(0.45, 1.2))
})


test_that("tabulation handles equal known components", {
  skip_if_not_installed("mockery")
  samples <- list(rnorm(20), rnorm(20))
  admixMod <- list(admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1)),
                   admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1)))
  fake_estim <- list(estimated_mixing_weights = c(0.4, 0.5), integ.supp = seq(-2, 2, length.out = 20), p.X.fixed = 0.3)

  mockery::stub(IBM_tabul_stochasticInteg, "estim_IBM", function(...) fake_estim)
  mockery::stub(IBM_tabul_stochasticInteg, "IBM_empirical_contrast", function(...) 0.5)
  mockery::stub(IBM_tabul_stochasticInteg, "detect_support_type", function(...) "Continuous")
  mockery::stub(IBM_tabul_stochasticInteg, "IBM_normalization_term", function(...) diag(3))
  mockery::stub(IBM_tabul_stochasticInteg, "estimVarCov_empProcess_Rcpp", function(...) diag(10))
  mockery::stub(IBM_tabul_stochasticInteg, "is_equal_knownComp", function(...) TRUE)
  mockery::stub(IBM_tabul_stochasticInteg, "sim_gaussianProcess", function(...) list(traj1 = rnorm(10)))

  res <- IBM_tabul_stochasticInteg(samples = samples, admixMod = admixMod, n_sim_tab = 2, parallel = FALSE)

  expect_true(is.numeric(res$U_sim))
})


test_that("tabulation is reproducible with seed", {
  set.seed(123)
  x1 <- rnorm(20)
  set.seed(123)
  x2 <- rnorm(20)
  expect_equal(x1, x2)
})


test_that("IBM_tabul_stochasticInteg couvre le support discret (multinom)", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.6, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.5, 0.3)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.4, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.5, 0.3)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1  <- admix_model(knownComp_dist  = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  mod2  <- admix_model(knownComp_dist  = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3)))
  res <- IBM_tabul_stochasticInteg(samples = list(data1, data2), admixMod = list(mod1, mod2), min_size = NULL,
                                   n.varCovMat = 10, n_sim_tab = 5, parallel = FALSE)

  expect_true(is.numeric(res$U_sim))
  expect_true(length(res$U_sim) > 0)
})
