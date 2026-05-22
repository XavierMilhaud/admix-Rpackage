test_that("multinom can only be mixed with multinom", {
  expect_error(
    twoComp_mixt(comp.dist = list("multinom", "norm"), comp.param = list(list(size = 1, prob = c(0.5,0.5)),
                                                                         list(mean = 0, sd = 1))),
    "multinom")
})

test_that("multinomial mixture returns discrete labels", {
  obj <- twoComp_mixt(n = 100, comp.dist = list("multinom", "multinom"),
                      comp.param = list(list(size = 1, prob = c(0.2,0.3,0.5)),
                                        list(size = 1, prob = c(0.4,0.4,0.2))))
  expect_true(all(obj$mixt.data %in% c(1,2,3)))
})

test_that("all output vectors have coherent sizes", {
  obj <- twoComp_mixt(n = 250, comp.dist = list("multinom", "multinom"),
                      comp.param = list(list(size = 1, prob = c(0.2,0.3,0.5)),list(size = 1, prob = c(0.5,0.3,0.2))))
  expect_length(obj$mixt.data, 250)
  expect_equal(length(obj$comp1.data) + length(obj$comp2.data),
               250)
})
