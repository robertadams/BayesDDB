test_that("bdb_single_arm_bin returns list with expected elements", {
  set.seed(1)
  n0 <- 40
  n1 <- 25
  X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
  Z <- c(rep(0, n0), rep(1, n1))
  Y <- c(rbinom(n0, 1, 0.4), rbinom(n1, 1, 0.35))
  out <- bdb_single_arm_bin(
    Y = Y,
    Z = Z,
    X = X,
    A = 10,
    s = 2,
    a = 1,
    K = 500,
    show_plots = FALSE
  )
  expect_type(out, "list")
  expect_true("theta_o" %in% names(out))
  expect_true("A_eff" %in% names(out))
  expect_true("Width" %in% names(out))
  expect_true("CrI" %in% names(out))
  expect_length(out$theta_o, 3)
})

test_that("bdb_single_arm_bin with double_adj = FALSE runs", {
  set.seed(2)
  n0 <- 30
  n1 <- 20
  X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
  Z <- c(rep(0, n0), rep(1, n1))
  Y <- c(rbinom(n0, 1, 0.3), rbinom(n1, 1, 0.35))
  out <- bdb_single_arm_bin(
    Y = Y,
    Z = Z,
    X = X,
    A = 5,
    s = 2,
    double_adj = FALSE,
    K = 300,
    show_plots = FALSE
  )
  expect_type(out, "list")
  expect_true("theta_o" %in% names(out))
})
