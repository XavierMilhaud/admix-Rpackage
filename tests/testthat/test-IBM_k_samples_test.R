test_that("IBM_k_samples_test fails with invalid admixMod", {
  fake_mod <- list(a = 1)
  expect_error(
    IBM_k_samples_test(samples = list(rnorm(10), rnorm(10), rnorm(10)), admixMod = list(fake_mod, fake_mod, fake_mod)),
    "admixMod")
})

test_that("icv testing requires at least two samples", {
  mod <- admix_model("norm", list(mean = 0, sd = 1))
  expect_error(
    admix_test(samples = list(rnorm(100)), admixMod = list(mod), test_method = "icv"),
    "requires at least TWO samples"
  )
})


test_that("admix_test dispatches to IBM_k_samples_test", {
  skip_on_cran()

  fake_test <- list(statistic = 3.4, p.value = 0.01)
  local_mocked_bindings(
    IBM_k_samples_test = function(samples, admixMod, conf_level, ...) { fake_test }
  )
  mod1 <- admix_model("norm", list(mean = 0, sd = 1))
  mod2 <- admix_model("exp", list(rate = 1))
  res <- admix_test(samples = list(rnorm(100), rexp(120)), admixMod = list(mod1, mod2), test_method = "icv")

  expect_s3_class(res, "IBM_test")
  expect_s3_class(res, "htest")
  expect_equal(res$statistic, 3.4)
})


test_that("IBM_k_samples_test fails with invalid admixMod", {
  skip_on_cran()

  obj <- create_test_samples()
  bad_admix <- list(1, 2)
  expect_error(
    IBM_k_samples_test(samples = obj$samples, admixMod = bad_admix),
    "Argument 'admixMod' is not correctly specified"
  )
})


test_that("IBM_k_samples_test dispatches to IBM_2samples_test for K=2", {
  skip_on_cran()

  obj <- create_test_samples()
  res <- IBM_k_samples_test(samples = obj$samples, admixMod = obj$admix, n_sim_tab = 5)

  expect_s3_class(res, "IBM_test")
  expect_equal(res$n_populations, 2)
})


test_that("IBM_k_samples_test works with K > 2", {
  skip_on_cran()

  set.seed(1)
  obj <- create_test_samples()
  samples3 <- c(obj$samples, list(rnorm(110)))
  admix3 <- c(obj$admix, list(admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))))
  res <- IBM_k_samples_test(samples = samples3, admixMod = admix3, tune_penalty = FALSE, n_sim_tab = 5, parallel = FALSE)

  expect_equal(res$n_populations, 3)
  expect_true(is.matrix(res$discrepancy_matrix))
  expect_true(is.numeric(res$p.value))
})


test_that("IBM_k_samples_test fonctionne avec parallel = TRUE (2 échantillons)", {
  skip_on_cran()
  skip_if_not_installed("doParallel")
  skip_if_not_installed("doRNG")

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 250, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 1,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 1, sd = 1))
  res <- suppressMessages(
    IBM_k_samples_test(samples = list(data1, data2), admixMod = list(mod1, mod2), conf_level = 0.95,
                       n_sim_tab = 5, parallel = TRUE, n_cpu = 2)
  )

  expect_s3_class(res, "IBM_test")
})


test_that("IBM_k_samples_test couvre tune_penalty = TRUE avec K = 3", {
  skip_on_cran()   # lent

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 380, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 1,  sd = 1)))
  mixt3 <- twoComp_mixt(n = 350, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 2,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  data3 <- get_mixture_data(mixt3)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 1, sd = 1))
  mod3  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 2, sd = 1))
  res <- suppressMessages(
    IBM_k_samples_test(samples = list(data1, data2, data3), admixMod = list(mod1, mod2, mod3), conf_level = 0.95,
                       n_sim_tab = 5, tune_penalty = TRUE, parallel = FALSE)
  )

  expect_s3_class(res, "IBM_test")
  expect_true(!is.na(res$tuning_param["Tuned Gamma"]))
  expect_true(!is.na(res$tuning_param["Tuned C"]))
  expect_true(is.logical(res$penalty_nullHyp))
})


test_that("IBM_k_samples_test utilise sim_U quand il est fourni", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 250, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 1,  sd = 1)))
  mixt3 <- twoComp_mixt(n = 280, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 2,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  data3 <- get_mixture_data(mixt3)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 1, sd = 1))
  mod3  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 2, sd = 1))
  # Fournir sim_U pré-calculé pour éviter la tabulation
  sim_U_precomputed <- abs(rnorm(200, mean = 0.5, sd = 0.3))
  res <- suppressMessages(
    IBM_k_samples_test(samples = list(data1, data2, data3), admixMod = list(mod1, mod2, mod3),
                       conf_level = 0.95, n_sim_tab = 5, tune_penalty = FALSE, sim_U = sim_U_precomputed)
  )

  expect_s3_class(res, "IBM_test")
  expect_equal(res$tabulated_dist, sim_U_precomputed)
})


test_that("IBM_k_samples_test stoppe si sim_U est entièrement NA", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 250, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 1,  sd = 1)))
  mixt3 <- twoComp_mixt(n = 280, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 2,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  data3 <- get_mixture_data(mixt3)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 1, sd = 1))
  mod3  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 2, sd = 1))

  expect_error(
    suppressMessages(
      IBM_k_samples_test(samples = list(data1, data2, data3), admixMod = list(mod1, mod2, mod3), conf_level = 0.95,
                         n_sim_tab = 5, tune_penalty = FALSE, sim_U = c(NA, NA, NA))   # tous NA
    ),
    "tabulate"
  )
})
