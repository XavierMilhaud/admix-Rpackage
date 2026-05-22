test_that("print.admix_model displays correct information", {
  obj <- admix_model("norm", list(mean = 0, sd = 1))
  expect_output(print(obj), "Known component distribution: norm")
  expect_output(print(obj), "mean = 0")
  expect_output(print(obj), "sd = 1")
})

test_that("print.admix_model returns object invisibly", {
  obj <- admix_model("norm", list(mean = 0, sd = 1))
  expect_invisible(print(obj))
})

test_that("summary.admix_model prints summary information", {
  obj <- admix_model("pois", list(lambda = 2))
  expect_output(summary(obj), "Known component distribution: pois")
})
