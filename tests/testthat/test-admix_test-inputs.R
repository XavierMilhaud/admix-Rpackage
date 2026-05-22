test_that("admix_test requires samples and admixMod to be lists", {
  expect_error(
    admix_test(samples = 1:10, admixMod = list()),
    "Please provide sample\\(s\\) AND admixture model\\(s\\) in a list"
  )
  expect_error(
    admix_test(samples = list(1:10), admixMod = 1),
    "Please provide sample\\(s\\) AND admixture model\\(s\\) in a list"
  )
})

test_that("admix_test validates admix_model objects", {
  expect_error(
    admix_test(samples = list(rnorm(10)), admixMod = list(list(a = 1))),
    "Argument 'admixMod' is not correctly specified"
  )
})
