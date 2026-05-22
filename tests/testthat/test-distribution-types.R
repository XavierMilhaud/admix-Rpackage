test_that("continuous distributions are correctly identified", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("norm", "exp"), comp.param = list(list(mean = 0, sd = 1),
                                                                                  list(rate = 1)))
  expect_equal(obj$dist.type, c("Continuous", "Continuous"))
})

test_that("discrete distributions are correctly identified", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("pois", "geom"), comp.param = list(list(lambda = 3),
                                                                                   list(prob = 0.4)))
  expect_equal(obj$dist.type, c("Discrete", "Discrete"))
})

test_that("multinomial distributions are identified as multivariate", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("multinom", "multinom"),
                      comp.param = list(list(size = 1, prob = c(0.2,0.3,0.5)),
                                        list(size = 1, prob = c(0.4,0.4,0.2))))
  expect_equal(obj$dist.type, c("Multivariate", "Multivariate"))
})
