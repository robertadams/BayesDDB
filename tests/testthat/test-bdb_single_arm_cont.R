test_that("bdb_single_arm_cont runs when rjags is available", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  set.seed(1)
  n0 <- 30
  n1 <- 20
  X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
  Z <- c(rep(0, n0), rep(1, n1))
  Y <- c(rnorm(n0, 0.5, 1), rnorm(n1, 0.4, 1))
  out <- bdb_single_arm_cont(
    Y = Y,
    Z = Z,
    X = X,
    A = 5,
    s = 2,
    a = 1,
    n_iter = 100,
    show_plots = FALSE
  )
  expect_type(out, "list")
  expect_true("theta_o" %in% names(out))
  expect_true("A_eff" %in% names(out))
})
