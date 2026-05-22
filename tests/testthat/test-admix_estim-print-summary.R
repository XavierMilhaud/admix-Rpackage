test_that("print.admix_estim works for PS objects", {
  obj <- structure(
    list(estim_objects = list(list(estimated_mixing_weights = 0.75, population_sizes = 100)),
         sample_names = "Sample_1",
         call = quote(admix_estim())
    ),
    class = c("estim_PS", "admix_estim")
  )

  expect_output(admix:::print.admix_estim(obj), "Method: PS")
  expect_output(admix:::print.admix_estim(obj), "0.750")
  expect_invisible(admix:::print.admix_estim(obj))
})


test_that("summary.admix_estim displays sample summaries", {
  fake_obj <- structure(
    list(estim_objects = list(structure(list(), class = "estim_PS")),
         sample_names = "Sample_A",
         call = quote(admix_estim())
    ),
    class = c("estim_PS", "admix_estim")
  )

  local_mocked_bindings(
    summary.estim_PS = function(object, ...) { cat("fake summary") }
  )

  expect_output(admix:::summary.admix_estim(fake_obj), "Sample: Sample_A")
  expect_output(admix:::summary.admix_estim(fake_obj), "fake summary")
})


test_that("print.admix_estim couvre la branche estim_IBM (composantes distinctes)", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 5, sd = 2)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 5, sd = 2))
  res <- suppressMessages(admix_estim(samples = list(data1, data2), admixMod = list(mod1, mod2), est_method = "IBM"))

  out <- capture.output(print.admix_estim(res))
  expect_true(any(grepl("IBM", out, ignore.case = TRUE)))
  expect_true(any(grepl("mix_weight", out, ignore.case = TRUE)))
})


test_that("print.admix_estim couvre la branche estim_IBM (composantes égales)", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 5, sd = 2), list(mean = 0, sd = 1)))  # même composante connue
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))  # identique à mod1
  res <- suppressMessages(admix_estim(samples = list(data1, data2), admixMod = list(mod1, mod2), est_method = "IBM"))

  out <- capture.output(print.admix_estim(res))
  # Dans ce cas, la colonne "mix_weight_1st (fixed)" doit apparaître
  expect_true(any(grepl("fixed", out, ignore.case = TRUE)))
})


test_that("print.admix_estim affiche les variances pour estim_BVdk", {
  skip_on_cran()
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- suppressMessages(admix_estim(samples = list(data1), admixMod = list(mod1), est_method = "BVdk", compute_var = TRUE))
  out <- capture.output(print.admix_estim(res))

  expect_true(any(grepl("var", out, ignore.case = TRUE)))
})


test_that("print.admix_estim gère sample_names NULL", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- suppressMessages(admix_estim(samples = list(data1), admixMod = list(mod1), est_method = "PS"))
  # Forcer sample_names à NULL pour couvrir la branche défensive
  res$sample_names <- NULL
  out <- capture.output(print.admix_estim(res))

  expect_true(any(grepl("Sample_1", out)))  # fallback "Sample_k" attendu
})


test_that("print.admix_estim gère sample_names avec chaîne vide", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- suppressMessages(admix_estim(samples = list(data1), admixMod = list(mod1), est_method = "PS"))
  res$sample_names <- ""   # chaîne vide
  out <- capture.output(admix:::print.admix_estim(res))

  expect_true(any(grepl("Sample_1", out)))
})


test_that("summary.admix_estim couvre la branche estim_IBM", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.5, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  mixt2 <- twoComp_mixt(n = 600, weight = 0.7, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 5, sd = 2)))
  data1 <- get_mixture_data(mixt1)
  data2 <- get_mixture_data(mixt2)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  mod2  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 5, sd = 2))
  res <- suppressMessages(admix_estim(samples = list(data1, data2), admixMod = list(mod1, mod2), est_method = "IBM"))
  out <- capture.output(admix:::summary.admix_estim(res))

  expect_true(any(grepl("Pair", out, ignore.case = TRUE)))
  expect_true(any(grepl("IBM", out, ignore.case = TRUE)))
})


test_that("summary.admix_estim gère sample_names NULL", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- admix_estim(samples = list(data1), admixMod = list(mod1), est_method = "PS")
  res$sample_names <- NULL
  out <- capture.output(admix:::summary.admix_estim(res))

  expect_true(any(grepl("Sample_1", out)))
})
