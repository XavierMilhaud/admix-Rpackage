test_that("admix_estim requires samples and admixMod to be lists", {
  expect_error(admix_estim(samples = 1:10, admixMod = list()),
               "Please provide sample\\(s\\) AND admixture model\\(s\\) in a list")
  expect_error(admix_estim(samples = list(1:10), admixMod = 1),
               "Please provide sample\\(s\\) AND admixture model\\(s\\) in a list")
})


test_that("admix_estim validates admix_model objects", {
  expect_error(admix_estim(samples = list(1:10), admixMod = list(list(a = 1))),
               "Argument 'admixMod' is not correctly specified")
})


test_that("admix_estim retrieves sample names from function call", {
  fake_estim <- list(estimated_mixing_weights = 0.7, population_sizes = 100)
  local_mocked_bindings(
    estim_PS = function(samples, admixMod, ...) fake_estim,
    detect_support_type = function(x) "Continuous"
  )

  x1 <- rnorm(100)
  x2 <- rnorm(200)
  mod1 <- admix_model("norm", list(mean = 0, sd = 1))
  mod2 <- admix_model("norm", list(mean = 1, sd = 2))
  res <- admix_estim(samples = list(x1, x2), admixMod = list(mod1, mod2), est_method = "PS")

  expect_equal(res$sample_names, c("x1", "x2"))
})


test_that("admix_estim uses default sample names when names are unavailable", {
  fake_estim <- list(estimated_mixing_weights = 0.7, population_sizes = 100)
  local_mocked_bindings(
    estim_PS = function(samples, admixMod, ...) fake_estim,
    detect_support_type = function(x) "Continuous"
  )

  samples <- list(rnorm(100), rnorm(200))
  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_estim(samples = samples, admixMod = list(mod, mod), est_method = "PS")

  expect_equal(res$sample_names, c("Sample_1", "Sample_2"))
})


test_that("admix_estim stoppe si IBM est demandé avec un seul échantillon", {
  set.seed(1)
  mixt1 <- twoComp_mixt(n = 300, weight = 0.6, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = 2, sd = 0.5), list(mean = 0, sd = 1)))
  data1 <- get_mixture_data(mixt1)
  mod1  <- admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1))

  expect_error(
    admix_estim(samples = list(data1), admixMod = list(mod1), est_method = "IBM"),
    "two samples")
})


