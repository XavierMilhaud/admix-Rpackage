test_that("returned object has correct class and fields", {
  set.seed(123)
  obj <- twoComp_mixt(n = 100, weight = 0.7)
  expect_s3_class(obj, "twoComp_mixt")
  expect_named(obj,
               c("n","mixt.data","dist.type","comp.dist","comp.param","mix.prop","comp1.data","comp2.data","call")
  )
})

test_that("generated sample has expected size", {
  obj <- twoComp_mixt(n = 250)
  expect_length(obj$mixt.data, 250)
  expect_equal(length(obj$comp1.data) + length(obj$comp2.data), 250)
})

test_that("simulation is reproducible with set.seed", {
  set.seed(123)
  x1 <- twoComp_mixt(n = 100)
  set.seed(123)
  x2 <- twoComp_mixt(n = 100)
  expect_identical(x1$mixt.data, x2$mixt.data)
})

test_that("component data are subsets of mixture data", {
  obj <- twoComp_mixt(n = 500)
  expect_equal(length(obj$comp1.data) + length(obj$comp2.data),
               length(obj$mixt.data))
})

test_that("simulation runs without warnings", {
  expect_no_warning(twoComp_mixt(n = 100))
})

#test_that("Right simulation of the data", {
#  expect_equal(round(mean(twoComp_mixt(n = 100000,
#                                       weight = 0.7,
#                                       comp.dist = list("norm","norm"),
#                                       comp.param = list(list("mean" = 2, "sd" = 0.2),
#                                                         list("mean" = 4, "sd" = 0.1)))$mixt.data),1), 2.6)
#})
