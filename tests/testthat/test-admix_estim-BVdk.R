test_that("BVdk rejects discrete distributions", {
  local_mocked_bindings(
    detect_support_type = function(x) "Discrete"
  )

  mod <- admix_model("pois", list(lambda = 2))

  expect_error(
    admix_estim(samples = list(rpois(100, 2)), admixMod = list(mod), est_method = "BVdk"),
    "'BVdk' estimation method is not suitable to discrete random variables"
  )
})


test_that("BVdk prints symmetry assumption message", {
  fake_estim <- list(estimated_mixing_weights = 0.5, population_sizes = 100, estimated_locations = 0,
                     mix_weight_variance = NA, location_variance = NA)

  local_mocked_bindings(
    estim_BVdk = function(samples, admixMod, ...) fake_estim,
    detect_support_type = function(x) "Continuous"
  )

  mod <- admix_model("norm", list(mean = 0, sd = 1))

  expect_message(
    admix_estim(samples = list(rnorm(100)), admixMod = list(mod), est_method = "BVdk"),
    "symmetric probability density function"
  )
})
