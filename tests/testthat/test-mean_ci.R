test_that("mean_ci returns mean and interval", {
  x <- c(1, 2, 3, 4, 5)
  out <- BayesDDB:::mean_ci(x, 0.95)
  expect_named(out, c("mean", "lower", "upper"))
  expect_equal(out[["mean"]], 3)
  expect_true(out[["lower"]] < out[["mean"]])
  expect_true(out[["upper"]] > out[["mean"]])
})

test_that("mean_ci respects conf_level", {
  x <- rnorm(100, mean = 0, sd = 1)
  out95 <- BayesDDB:::mean_ci(x, 0.95)
  out99 <- BayesDDB:::mean_ci(x, 0.99)
  expect_true(out99[["lower"]] < out95[["lower"]])
  expect_true(out99[["upper"]] > out95[["upper"]])
})
