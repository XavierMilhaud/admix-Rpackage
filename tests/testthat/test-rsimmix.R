test_that("Right simulation of the data", {
  expect_equal(round(mean(rsimmix(n = 100000,
                                  unknownComp_weight = 0.7,
                                  comp.dist = list(f = 'norm', g = 'norm'),
                                  comp.param = list(f = list(mean = 2, sd = 0.2),
                                                    g = list(mean = 4, sd = 0.1)))[['mixt.data']]),1), 2.6)
})