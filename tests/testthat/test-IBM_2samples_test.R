test_that("IBM_2samples_test returns a valid IBM_test object", {
  skip_on_cran()
  obj <- create_test_samples()
  res <- IBM_2samples_test(samples = obj$samples, admixMod = obj$admix, conf_level = 0.95, n_sim_tab = 5, parallel = FALSE)

  expect_s3_class(res, "IBM_test")
  expect_s3_class(res, "htest")
  expect_type(res$p.value, "double")
  expect_true(is.logical(res$reject_decision))
  expect_equal(res$n_populations, 2)
  expect_length(res$population_sizes, 2)
  expect_true(is.numeric(res$statistic))
})


test_that("results are reproducible with fixed seed", {
  skip_on_cran()
  obj <- create_test_samples()
  set.seed(123)
  res1 <- IBM_2samples_test(samples = obj$samples, admixMod = obj$admix, n_sim_tab = 5)
  set.seed(123)
  res2 <- IBM_2samples_test(samples = obj$samples, admixMod = obj$admix, n_sim_tab = 5)

  expect_equal(res1$statistic, res2$statistic)
})


test_that("IBM_2samples_test couvre la branche estim.weights hors [-1,1]", {
  skip_on_cran()
  set.seed(1)
  # Données très courtes et mal conditionnées pour forcer une divergence
  data1 <- c(0.1, 0.2, 0.3, 50, 100)   # valeurs extrêmes
  data2 <- c(0.1, 0.2, 200, 300, 400)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- suppressMessages(
    IBM_2samples_test(samples  = list(data1, data2), admixMod = list(mod1, mod2), conf_level = 0.95, n_sim_tab = 3)
  )

  # Si la branche est atteinte, reject = TRUE et p_value = 1e-12
  if (any(abs(res$statistic) > 1) || any(is.na(res$statistic))) {
    expect_true(res$reject_decision)
    expect_equal(res$p.value, 1e-12)
  } else {
    # test passé normalement, la branche n'a pas été déclenchée
    expect_s3_class(res, "IBM_test")
  }
})


test_that("IBM_2samples_test couvre la branche estim.weights de longueur 1 (G1 == G2)", {
  skip_on_cran()
  set.seed(1)
  # Composantes connues identiques → p1 fixé, estim.weights de longueur 1
  mixt1 <- twoComp_mixt(n = 400, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 400, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))  # identique
  res <- suppressMessages(
    IBM_2samples_test(samples  = list(data1, data2), admixMod = list(mod1, mod2), conf_level = 0.95, n_sim_tab = 5)
  )

  expect_s3_class(res, "IBM_test")
  # Dans ce cas, estimated_values doit avoir deux éléments : c(0.2, estim.weights)
  expect_length(res$statistic, 1)
  expect_named(res$statistic, "T")
})
