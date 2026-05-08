test_that("bdb_two_arms_bin returns list with expected elements", {
  set.seed(1)
  n0 <- 40
  n1 <- 30
  X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
  Z <- c(rep(0, n0), rep(1, n1))
  TRT <- c(rep(0, n0), rbinom(n1, 1, 0.5))
  Y <- c(rbinom(n0, 1, 0.4), rbinom(n1, 1, 0.35))
  out <- bdb_two_arms_bin(Y = Y, Z = Z, TRT = TRT, X = X, A = 10, s = 2, a = 1, K = 500, show_plots = FALSE)
  expect_type(out, "list")
  expect_true("theta_o" %in% names(out))
  expect_true("mean.theta" %in% names(out))
  expect_true("A_eff" %in% names(out))
  expect_true("Width" %in% names(out))
})

test_that("bdb_two_arms_bin with double_adj = FALSE runs", {
  set.seed(2)
  n0 <- 30
  n1 <- 24
  X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
  Z <- c(rep(0, n0), rep(1, n1))
  TRT <- c(rep(0, n0), rbinom(n1, 1, 0.5))
  Y <- c(rbinom(n0, 1, 0.35), rbinom(n1, 1, 0.4))
  out <- bdb_two_arms_bin(Y = Y, Z = Z, TRT = TRT, X = X, A = 5, s = 2, double_adj = FALSE, K = 300, show_plots = FALSE)
  expect_type(out, "list")
  expect_true("theta_o" %in% names(out))
})
