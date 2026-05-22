test_that("n must be a positive scalar", {
  expect_error(twoComp_mixt(n = -1), "`n` must be a positive integer.")
  expect_error(twoComp_mixt(n = c(1,2)), "`n` must be a positive integer.")
  expect_error(twoComp_mixt(n = "10"), "`n` must be a positive integer.")
})

test_that("weight must belong to (0,1)", {
  expect_error(twoComp_mixt(weight = 0), "`weight` must belong to \\(0,1\\).")
  expect_error(twoComp_mixt(weight = 1), "`weight` must belong to \\(0,1\\).")
  expect_error(twoComp_mixt(weight = -0.2), "`weight` must belong to \\(0,1\\).")
})

test_that("exactly two component distributions are required", {
  expect_error(twoComp_mixt(comp.dist = list("norm"), comp.param = list(list(mean = 0, sd = 1))),
               "Please provide exactly two component distributions.")
})

test_that("n must be integer", {
  expect_error(twoComp_mixt(n = 1.5), "integer")
})
