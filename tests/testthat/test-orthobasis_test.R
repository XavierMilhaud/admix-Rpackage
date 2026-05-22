test_that("orthobasis_test stoppe si admixMod est mal spécifié", {
  set.seed(1)
  data1 <- rnorm(200)
  data2 <- rnorm(200)

  expect_error(
    orthobasis_test(samples = list(data1, data2), admixMod = list("not_a_model", "not_a_model"), support = "Real"),
    "admixMod")
})

test_that("orthobasis_test stoppe si nb_echBoot <= 1 avec méthode PS", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 200, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 200, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  expect_error(
    orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), est_method = "PS",
                    nb_echBoot = 1, support = "Real"),
    "bootstrap"
  )
})

test_that("orthobasis_test stoppe si s est hors ]0, 0.5[", {
  set.seed(1)
  data1 <- rnorm(200)
  data2 <- rnorm(200)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  expect_error(
    orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), s = 0, support = "Real"),
    "penalty"
  )
  expect_error(
    orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), s = 0.5, support = "Real"),
    "penalty"
  )
})


test_that("orthobasis_test couvre ask_poly_param = TRUE", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  # Mock readline pour simuler la saisie utilisateur
  local_mocked_bindings(
    readline = function(prompt = "") {
      if (grepl("K", prompt)) "3" else "0.25"
    },
    .package = "base"
  )
  res <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), ask_poly_param = TRUE, support = "Real")

  expect_s3_class(res, "orthobasis_test")
})


test_that("orthobasis_test fonctionne avec la méthode PS", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 1, sd = 1), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  res <- suppressMessages(
    orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2),
                    est_method = "PS", nb_echBoot = 10, support = "Real"))

  expect_s3_class(res, "orthobasis_test")
  expect_true(is.numeric(res$statistic))
  expect_true(is.numeric(res$p.value))
})


test_that("orthobasis_test rejette H0 quand les composantes inconnues sont clairement différentes", {
  skip_on_cran()

  set.seed(1)
  # Composantes inconnues très différentes : N(0,1) vs N(10,1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0,  sd = 1), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 10, sd = 1), list(mean = 0,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res_H1 <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), support = "Real")

  expect_true(res_H1$reject_decision)        # rej <- TRUE atteint
  expect_lt(res_H1$p.value, 0.05)
})


test_that("print.orthobasis_test s'exécute sans erreur", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0,  sd = 1), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0, sd = 1), list(mean = 0,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res_H0 <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), support = "Real")

  out <- capture.output(print.orthobasis_test(res_H0))
  expect_true(any(grepl("null hypothesis", out, ignore.case = TRUE)))
})


test_that("print.orthobasis_test affiche 'Yes' quand H0 est rejetée", {
  skip_on_cran()

  set.seed(1)
  # Composantes inconnues très différentes : N(0,1) vs N(10,1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0,  sd = 1), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 10, sd = 1), list(mean = 0,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res_H1 <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), support = "Real")
  out <- capture.output(print.orthobasis_test(res_H1))   # res_H1 = résultat du groupe D
  expect_true(any(grepl("Yes", out)))
})


test_that("print.orthobasis_test gère p-value arrondie à 0", {
  skip_on_cran()

  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0,  sd = 1), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 10, sd = 1), list(mean = 0,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res_H1 <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), support = "Real")
  # Forcer un objet avec p.value très petite
  res_tiny_p <- res_H1
  res_tiny_p$p.value <- 1e-15
  out <- capture.output(print.orthobasis_test(res_tiny_p))
  expect_true(any(grepl("1e-12", out)))
})


test_that("summary.orthobasis_test s'exécute sans erreur", {
  skip_on_cran()

  set.seed(1)
  # Composantes inconnues très différentes : N(0,1) vs N(10,1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0,  sd = 1), list(mean = 0,  sd = 1)))
  mixt2 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 0, sd = 1), list(mean = 0,  sd = 1)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2 <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res_H0 <- orthobasis_test(samples = list(data1, data2), admixMod = list(mod1, mod2), support = "Real")
  out <- capture.output(summary.orthobasis_test(res_H0))
  expect_true(any(grepl("samples", out, ignore.case = TRUE)))
  expect_true(any(grepl("statistic", out, ignore.case = TRUE)))
})
