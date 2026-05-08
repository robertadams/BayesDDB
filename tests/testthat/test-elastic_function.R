test_that("elastic_function returns numeric vector of same length as ppp", {
  ppp <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  out <- elastic_function(ppp, a = 1)
  expect_type(out, "double")
  expect_length(out, length(ppp))
})

test_that("elastic_function gives 1 at ppp = 0.5", {
  expect_equal(elastic_function(0.5, a = 1), 1)
  expect_equal(elastic_function(0.5, a = 5), 1)
})

test_that("elastic_function output is in [0, 1]", {
  x <- seq(0, 1, by = 0.1)
  out <- elastic_function(x, a = 1)
  expect_true(all(out >= 0 & out <= 1))
})

test_that("elastic_function errors on invalid ppp", {
  expect_error(elastic_function(-0.1, a = 1))
  expect_error(elastic_function(1.5, a = 1))
  expect_error(elastic_function("x", a = 1))
})

test_that("elastic_function errors on a < 1", {
  expect_error(elastic_function(0.5, a = 0.5))
})
