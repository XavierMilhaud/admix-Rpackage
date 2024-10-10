test_that("Right simulation of the data", {
  expect_equal(round(mean(twoComp_mixt(n = 100000,
                                       weight = 0.7,
                                       comp.dist = list("norm","norm"),
                                       comp.param = list(list("mean" = 2, "sd" = 0.2),
                                                         list("mean" = 4, "sd" = 0.1)))$mixt.data),1), 2.6)
})
