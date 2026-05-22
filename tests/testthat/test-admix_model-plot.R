test_that("plot.admix_model works for continuous distributions", {
  obj <- admix_model("norm", list(mean = 0, sd = 1))
  expect_no_error(plot(obj))
})

test_that("plot.admix_model works for discrete distributions", {
  obj <- admix_model("pois", list(lambda = 3))
  expect_no_error(plot(obj))
})

test_that("plot.admix_model works for multinomial distributions", {
  obj <- admix_model("multinom", list(size = 1, prob = c(0.2, 0.5, 0.3)))
  expect_no_error(plot(obj))
})
