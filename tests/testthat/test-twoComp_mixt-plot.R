test_that("plot method runs without error", {
  obj <- twoComp_mixt(n = 100)
  expect_no_error(plot(obj))
})

test_that("continuous plot branch works", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("norm", "exp"), comp.param = list(list(mean = 0, sd = 1),
                                                                                  list(rate = 1)))
  expect_no_error(plot(obj))
})

test_that("discrete plot branch works", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("pois", "geom"), comp.param = list(list(lambda = 3),
                                                                                   list(prob = 0.5)))
  expect_no_error(plot(obj))
})

test_that("multinomial plot branch works", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("multinom", "multinom"),
                      comp.param = list(list(size = 1, prob = c(0.3,0.3,0.4)),list(size = 1, prob = c(0.5,0.2,0.3))))
  expect_no_error(plot(obj))
})

test_that("add_plot = TRUE works correctly", {
  obj1 <- twoComp_mixt(n = 100)
  obj2 <- twoComp_mixt(n = 100, comp.dist = list("exp", "norm"), comp.param = list(list(rate = 1),list(mean = 2, sd = 1)))
  plot(obj1)
  expect_no_error(plot(obj2, add_plot = TRUE))
})
