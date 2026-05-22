test_that("gaussianity_test fails with invalid admixMod", {
  expect_error(
    gaussianity_test(sample = rnorm(20), admixMod = list()),
    "Argument 'admixMod' is not correctly specified")
})


test_that("interactive parameter branch works", {
  vals <- c("3", "0.25")
  i <- 0

  local_mocked_bindings(
    readline = function(...) {
      i <<- i + 1
      vals[i]
    },
    .package = "base"
  )

  samples <- rnorm(30)
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- gaussianity_test(sample = samples, admixMod = admixMod, ask_poly_param = TRUE)

  expect_s3_class(res, "gaussianity_test")
})


test_that("multinomial contamination model is rejected", {
  admixMod <- structure(
    list(comp.dist = list(known = "multinom"), comp.param = list(known = list(size = 1))),
    class = "admix_model"
  )

  expect_error(
    gaussianity_test(sample = rnorm(20), admixMod = admixMod),
    "multinomial known distribution is not supported")
})


test_that("invalid s parameter is rejected", {
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))

  expect_error(
    gaussianity_test(sample = rnorm(20), admixMod = admixMod, s = 0),
    "penalty exponent"
  )
  expect_error(
    gaussianity_test(sample = rnorm(20), admixMod = admixMod, s = 0.5),
    "penalty exponent"
  )
})


test_that("cubature fallback branch is covered", {
  skip_if_not_installed("cubature")

  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  local_mocked_bindings(
    integrate = function(...) stop("boom"),
    .package = "stats"
  )
  res <- gaussianity_test(sample = rnorm(30), admixMod = admixMod)

  expect_s3_class(res, "gaussianity_test")
})


test_that("reject branch is covered", {
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  local_mocked_bindings(
    qchisq = function(...) -Inf,
    .package = "stats"
  )
  res <- gaussianity_test(sample = rnorm(50), admixMod = admixMod)

  expect_true(res$reject_decision)
})


test_that("print.gaussianity_test prints output", {
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- gaussianity_test(sample = rnorm(40), admixMod = admixMod)

  expect_output(print.gaussianity_test(res), "p-value")
})


test_that("print.gaussianity_test handles tiny p-values", {
  obj <- list(p.value = 0, reject_decision = TRUE, call = quote(test()))
  class(obj) <- "gaussianity_test"

  expect_output(print.gaussianity_test(obj), "1e-12")
})


test_that("summary.gaussianity_test prints summary", {
  admixMod <- admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1))
  res <- gaussianity_test(sample = rnorm(40), admixMod = admixMod)

  expect_output(summary.gaussianity_test(res), "Test statistic")
  expect_output(summary.gaussianity_test(res), "Estimates")
})


test_that("summary handles zero p-values", {
  obj <- list(
    call = quote(test()),
    population_size = 10,
    admixture_model = admix_model(knownComp_dist = "norm", knownComp_param = list(mean = 0, sd = 1)),
    method = "Gaussianity test",
    reject_decision = TRUE,
    confidence_level = 0.95,
    p.value = 0,
    selected_rank = 1,
    statistic = 100,
    var_statistic = diag(3),
    estimate = c(0.5, 0, 1)
  )
  class(obj) <- c("gaussianity_test", "htest")

  expect_output(summary.gaussianity_test(obj), "1e-12")
})
