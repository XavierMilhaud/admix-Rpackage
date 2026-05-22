test_that("admix_estim returns PS estimator object", {
  fake_estim <- list(estimated_mixing_weights = 0.7, population_sizes = 100)

  local_mocked_bindings(
    estim_PS = function(samples, admixMod, ...) fake_estim,
    detect_support_type = function(x) "Continuous"
  )

  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_estim(samples = list(rnorm(100)), admixMod = list(mod), est_method = "PS")

  expect_s3_class(res, "estim_PS")
  expect_s3_class(res, "admix_estim")
  expect_length(res$estim_objects, 1)
  expect_equal(res$estim_objects[[1]]$estimated_mixing_weights, 0.7)
})


test_that("PS estimation works on simple simulated data", {
  set.seed(123)
  x <- c(rnorm(200, -2, 1), rnorm(100, 0, 1))
  mod <- admix_model("norm", list(mean = 0, sd = 1))
  res <- admix_estim(samples = list(x), admixMod = list(mod), est_method = "PS")
  expect_s3_class(res, "estim_PS")
  w <- res$estim_objects[[1]]$estimated_mixing_weights
  expect_true(w > 0)
  expect_true(w < 1)
})
