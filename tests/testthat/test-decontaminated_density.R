# ==============================================================================
# Helpers locaux
# ==============================================================================

make_decontam_continuous <- function() {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 400, weight = 0.4, comp.dist  = list("norm", "norm"),
                        comp.param = list(list(mean = -2, sd = 0.5), list(mean = 0,  sd = 1)))
  list(data = get_mixture_data(mixt1),
       admixMod = admix_model(knownComp_dist  = "norm", knownComp_param = list(mean = 0, sd = 1)))
}

make_decontam_discrete <- function() {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 1000, weight = 0.6, comp.dist  = list("pois", "pois"),
                        comp.param = list(list(lambda = 3), list(lambda = 2)))
  list(data = get_mixture_data(mixt1),
       admixMod = admix_model(knownComp_dist  = "pois", knownComp_param = list(lambda = 2)))
}

make_decontam_multinom <- function() {
  set.seed(42)
  mixt1 <- twoComp_mixt(n = 500, weight = 0.7, comp.dist  = list("multinom", "multinom"),
                        comp.param = list(list(size = 1, prob = c(0.2, 0.3, 0.5)),
                                          list(size = 1, prob = c(0.1, 0.6, 0.3))))
  list(data = get_mixture_data(mixt1),
       admixMod = admix_model(knownComp_dist  = "multinom", knownComp_param = list(size = 1, prob = c(0.1, 0.6, 0.3))))
}


test_that("decontaminated_density returns correct object for continuous data", {
  set.seed(123)
  sample1 <- rnorm(200)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  obj <- decontaminated_density(sample1 = sample1, admixMod = admixMod, estim.p = 0.8)

  expect_s3_class(obj, "decontaminated_density")
  expect_equal(obj$support, "Continuous")
  expect_true(is.function(obj$decontaminated_density_fun))

  vals <- obj$decontaminated_density_fun(c(-1, 0, 1))
  expect_type(vals, "double")
  expect_length(vals, 3)
  expect_true(all(vals >= 0))
})


test_that("decontaminated_density works for discrete data", {
  set.seed(123)
  sample1 <- rpois(300, lambda = 2)
  admixMod <- admix_model(knownComp_dist = "pois", knownComp_param = list(lambda = 1))
  obj <- decontaminated_density(sample1 = sample1, admixMod = admixMod, estim.p = 0.7)

  expect_equal(obj$support, "Discrete")

  vals <- obj$decontaminated_density_fun(c(0,1,2))
  expect_length(vals, 3)
  expect_true(all(vals >= 0))
})


test_that("decontaminated_density validates inputs", {
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  expect_error(
    decontaminated_density(sample1 = "bad", admixMod = admixMod, estim.p = 0.5),
    "'sample1' and 'estim.p' must be numeric")

  expect_error(
    decontaminated_density(sample1 = rnorm(10), admixMod = "bad", estim.p = 0.5),
    "must be of class 'admix_model'")

  expect_error(
    decontaminated_density(sample1 = rnorm(10), admixMod = admixMod, estim.p = 0),
    "'estim\\.p' must be in \\(0,1\\]\\.")

  expect_error(
    decontaminated_density(sample1 = rnorm(10), admixMod = admixMod, estim.p = 1.5),
    "'estim\\.p' must be in \\(0,1\\]\\.")
})


test_that("decontaminated_density : branche multinom/approxfun couverte", {
  d <- make_decontam_multinom()
  res <- decontaminated_density(sample1  = d$data, admixMod = d$admixMod, estim.p  = 0.7)
  expect_s3_class(res, "decontaminated_density")

  vals <- res$decontaminated_density_fun(1:3)
  expect_true(is.numeric(vals))
  expect_true(all(vals >= 0))
})


test_that("summary.decontaminated_density : branche support discret", {
  d <- make_decontam_discrete()
  res <- decontaminated_density(sample1  = d$data, admixMod = d$admixMod, estim.p  = 0.6)

  expect_equal(res$support, "Discrete")
  out <- capture.output(summary(res))
  expect_true(any(grepl("Count table", out, ignore.case = TRUE)))
})


test_that("plot.decontaminated_density : branche x_val non NULL couverte", {
  d <- make_decontam_continuous()
  res <- decontaminated_density(sample1  = d$data, admixMod = d$admixMod, estim.p  = 0.4)
  x_val_custom <- seq(-3, 3, length.out = 50)

  expect_silent(plot(res, x_val = x_val_custom))
})


test_that("decontaminated_cdf : branche pmultinom -> stepfun couverte", {
  d <- make_decontam_multinom()
  F1 <- decontaminated_cdf(sample1  = d$data, admixMod = d$admixMod, estim.p  = 0.7)
  expect_true(is.function(F1))

  vals <- F1(1:3)
  expect_true(is.numeric(vals))
  expect_true(all(vals >= 0))
  expect_true(all(vals <= 1))
  expect_true(all(diff(vals) >= 0))   # croissante
})
