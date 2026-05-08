
## Not run: 
# Simulate example data
library(patchwork)

set.seed(123)
N0 <- 80   # external subjects
N1 <- 40   # current trial (single arm)

# Covariates
X <- rbind(matrix(rnorm(N0 * 2), ncol = 2), matrix(rnorm(N1 * 2), ncol = 2))
Z <- c(rep(0, N0), rep(1, N1))
Y <- c(rnorm(N0, 50, 5), rnorm(N1, 52, 5))

# Run BDB analysis (requires JAGS)
result <- bdb_single_arm_cont(Y = Y, Z = Z, X = X,
                              A = 20, s = 3, a = 1,
                              za = 0.05, n_iter = 1000,
                              double_adj = TRUE, show_plots = TRUE)

result$plots$rates_by_stratum +
  result$plots$overall_rates +
  result$plots$mcmc +
  result$plots$ps_stratum +
  result$plots$ps_group +
  result$plots$forest
  
  
  
  

  