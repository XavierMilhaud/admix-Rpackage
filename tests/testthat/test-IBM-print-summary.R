test_that("print method works", {
  skip_on_cran()

  obj <- create_test_samples()
  res <- IBM_2samples_test(samples = obj$samples, admixMod = obj$admix, n_sim_tab = 5)

  expect_output(print(res), "p-value")
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


test_that("print.IBM_test affiche '1e-12' quand p-value arrondie à 0", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 350, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 1,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 1, sd = 1))
  res <- suppressMessages(
    IBM_2samples_test(samples  = list(data1, data2), admixMod = list(mod1, mod2), conf_level = 0.95, n_sim_tab = 5)
  )
  # Forcer p.value très petite pour déclencher la branche
  res$p.value <- 1e-15
  out <- capture.output(print.IBM_test(res))

  expect_true(any(grepl("1e-12", out)))
})

test_that("summary.IBM_test s'exécute sans erreur", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 350, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 3, sd = 0.5), list(mean = 1,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 1, sd = 1))
  res <- suppressMessages(
    IBM_2samples_test(samples  = list(data1, data2), admixMod = list(mod1, mod2), conf_level = 0.95, n_sim_tab = 5)
  )
  out <- capture.output(summary.IBM_test(res))

  expect_true(any(grepl("samples", out, ignore.case = TRUE)))
  expect_true(any(grepl("statistic", out, ignore.case = TRUE)))
  expect_true(any(grepl("p-value", out, ignore.case = TRUE)))
})
