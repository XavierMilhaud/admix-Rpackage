test_that("print method works", {
  obj <- twoComp_mixt(n = 50)
  expect_output(print(obj), "Number of observations")
})

test_that("print uses S3 dispatch", {
  obj <- twoComp_mixt(n = 10)
  expect_true(methods::is(obj, "twoComp_mixt") || inherits(obj, "twoComp_mixt"))
})

test_that("summary method returns output", {
  obj <- twoComp_mixt(n = 10)
  expect_output(summary(obj), "Statistics related")
})
