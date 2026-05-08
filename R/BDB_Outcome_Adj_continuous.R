#' Outcome-adjustment-only helper for continuous two-arm case (internal)
#'
#' Applies outcome adjustment for Bayesian Dynamic Borrowing in the two-arm continuous setting.
#' Used internally; not called by the main \code{bdb_*} API. Requires JAGS.
#'
#' @param Y Numeric vector of continuous outcomes.
#' @param Z Numeric vector: 1 = current trial, 0 = external.
#' @param TRT Numeric vector: 1 = treatment, 0 = control.
#' @param X Matrix or data frame of baseline covariates.
#' @param A Numeric. Intended number of subjects to borrow.
#' @param a Numeric. Tuning parameter for elastic function.
#' @param za Numeric. Significance level for credible intervals.
#' @param n.iter Integer. Number of MCMC iterations.
#' @return List with estimate, credible interval, and diagnostic quantities.
#' @keywords internal
#' @noRd
BDB_Outcome_Adj_continuous <- function(Y, Z, TRT, X, A, a, za, n.iter) {
  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' is required. Install it and ensure JAGS is installed.", call. = FALSE)
  for (pkg in c("dplyr", "overlapping", "ggplot2", "ggpubr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Please install it.", call. = FALSE)
  }
  # BDB_binary() applies Bayesian Dynamic Borrowing method for binary outcomes 
  # 
  # References: 
  #
  # Farjat, A., Ji, Y., Kaiser, A., Lu, C., Pap, A., Potts, J., Wang, M. 
  # A Comprehensive Bayesian Double-Adjustment Approach to Dynamic Borrowing of 
  # External Data. Manuscript in progress.
  #
  # Wang, C., Li, H., Chen, W.C., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2019. 
  # Propensity score-integrated power prior approach for incorporating real-world 
  # evidence in single-arm clinical studies. 
  # Journal of biopharmaceutical statistics, 29(5), pp.731-748.
  #  
  # INPUT
  # Y: n-by-1 vector of binary outcomes (0s and 1s) 
  # Z: n-by-1 external data indicator (current: 1, external:0)
  # TRT: n-by-1 active treatment indicator (active treatment: 1, control: 0) 
  # X: n-by-p matrix of confounders
  # A: number of subjects intended to be borrowed from external data source
  # a: tuning parameter for arctangent elastic function
  # za: confidence level value
  # n.iter: number of samples from posterior distribution
  #
  # OUTPUT
  # n_1: number of current subjects  
  # n_0: number of external subjects
  # n_11: current subjects in treatment arm 
  # n_10: number of internal control subjects 
  # R_ITA: overall mean and 95% CrI of internal treatment arm
  # R_ICA: overall mean and 95% CrI of internal control arm
  # R_ECA: overall mean and 95% CrI of external control arm
  # A_eff: Overall effective number of subjects borrowed from external sources
  # PPP: posterior predictive probability
  # gPPP: transformed posterior predictive probability used for further discounting
  # gamma: power discounting parameter
  # LCrI: length of credible interval from point estimate
  # Estimate: 1-by-3 vector of treatment effect and corresponding 95% CrI. 
  # The treatment effect is defined as the difference of means between treatment 
  # and control arms 
  #
  # Plots
  # Figure 1: Trace plot, density and ACF from posterior samples of of rate difference
  # Figure 2: Overall mean and 95% CI for internal (ICA) and external control arm (ECA)
  ######################################################################################################
  # Function begins
  # Define data object
  data<-data.frame(Y, Z, TRT, X)
  
  ##################################################################################################
  ##################################################################################################
  # Bayesian model for the calculation of predictive probability of observing current mean given external data
  
  model_string_PPP<-"
model
{
  for (i in 1:n){
    y[i] ~ dnorm(m,pre)
  }
  m ~ dnorm(0,0.00001)
  pre ~ dgamma(0.01,0.01)
  new ~ dnorm(m,pre*N1)
}
"
  
  
  model_string_TREATMENT<-"
model
{
  for (i in 1:n){
    y[i] ~ dnorm(mu,pre)
  }
  mu ~ dnorm(0,0.00001)
  pre ~ dgamma(0.01,0.01)
}
"
  
  modelfile <- tempfile(fileext = ".txt")
  on.exit(unlink(modelfile), add = TRUE)
  cat(
    "
model{
  for(i in 1:n){
    phi[i] <- -(log(pre)/2-(y[i]-mu)^2*pre/2)*gamma[i] + 10000 # stabilizing constant
    zeros[i] ~ dpois(phi[i])
  }

# Prior distributions
  mu ~ dnorm(0, 0.00001)
  pre ~ dgamma(0.01, 0.01)
}
", file = modelfile)

  ######################################################################################################
  ######################################################################################################
  # Step: Organize data 
  ######################################################################################################
  ######################################################################################################
  # Sample size external data source
  n_0<-sum(data$Z==0)
  # Sample size current trial
  n_1<-sum(data$Z==1)
  # Sample size current treatment arm
  n_11<-sum(data$Z==1 & data$TRT==1)
  # Sample size current control
  n_10<-sum(data$Z==1 & data$TRT==0)
  # outcome for current treatment 
  y_11<-data$Y[data$Z==1 & data$TRT==1]
  # outcome for internal control
  y_10<-data$Y[data$Z==1 & data$TRT==0] 
  # Number of events for external data
  y_0<-data$Y[data$Z==0] 
  
  # Mean and 95% CI for ITA, ICA, ECA
  R_ITA<-mean_CI(y_11, 0.05)
  R_ICA<-mean_CI(y_10, 0.05)
  R_ECA<-mean_CI(y_0, 0.05)
  
  ######################################################################################################
  # Step: Calibration
  ######################################################################################################
  ######################################################################################################
  # Estimate predictive posterior probability by stratum
  #observed mean in internal control 
  e.m<-mean(y_10) 
  
  # JAGS model  
  model_PPP<-textConnection(model_string_PPP)
  # Initial values of parameters
  inits<-list(m=e.m)
  # Data list as input to JAGS
  #n_0<-sum(data$Z==0)
  #y_0<-data$Y[data$Z==0]
  #n_10<-sum(data$Z==1 & data$TRT==0)
  dataList<-list("n"=n_0,"y"=y_0,"N1"=n_10) 
    
  # Run MCMC by JAGS (sampling from external data assuming precision from current data)
  jagsfit <- rjags::jags.model(model_PPP, data=dataList, inits=inits, n.chains=1, n.adapt=1000) 
  # Draw MCMC samples
  params<-c("m","new") # Parameters to be monitored in MCMC
  out <- rjags::coda.samples(jagsfit, params, n.iter=n.iter)
  out<-do.call(rbind.data.frame, out)
  # Calculate posterior predictive probability
  PPP<-mean(out$new<e.m)
  
  ######################################################################################################
  ######################################################################################################
  # Step: Mapping 
  ######################################################################################################
  ######################################################################################################
  # a: value of the tuning parameter in the arctangent elastic function
  # Further discounting by the elastic function 
  gPPP<-atan(a*sin(PPP*pi))/atan(a) # mapped predictive posterior probability
  
  ######################################################################################################
  ######################################################################################################
  # Step 9: Discounting 
  ######################################################################################################
  ######################################################################################################
  # A: number of subjects intended to be borrowed from external source
  # Define discounting factor
  gamma<-min(1,A/n_0)
  
  ######################################################################################################
  ######################################################################################################
  # Outcome adjustment 
  ######################################################################################################
  ######################################################################################################
  # Treatment arm
  # JAGS model  
  model_TREAT<-textConnection(model_string_TREATMENT)
  # Initial values
  mu_ini<-mean(y_11)
  pre_ini<-1/var(y_11)
  inits<-list(mu=mu_ini, pre=pre_ini)
  # Data list as input to JAGS
  dataList<-list("n"=n_11,"y"=y_11) 
  # Run MCMC by JAGS
  fit_treat <- rjags::jags.model(model_TREAT, data=dataList, inits=inits, n.chains=1, n.adapt=10000) 
  # Draw MCMC samples
  params_TREAT<-c("mu","pre") # Parameters to be monitored in MCMC
  out <- rjags::coda.samples(fit_treat, variable.names = params_TREAT, n.iter=n.iter)
  out<-do.call(rbind.data.frame, out)
  # Posterior samples
  mu_TREAT<-out$mu
  
  # Internal and external controls
  # Data list
  n_<-sum(data$TRT==0)
  y_<-data$Y[data$TRT==0]
  gamma_<-data$Z + (data$Z==0)*gamma*gPPP
  zeros_<-rep(0,n_) 
  data_CONT<-list(n=n_, y=y_, gamma=gamma_, zeros=zeros_)
  # Initial values
  mu_<-mean(y_)
  pre_<-1/var(y_)
  inits_CONT<-function(){ list( mu=mu_, pre=pre_) }
  # Fit model 
  fit_cont <- rjags::jags.model(file = modelfile, data = data_CONT, inits = inits_CONT, n.chains = 1)
  update(fit_cont, n.iter)
  # Parameters to monitor during MCMC
  params_CONT<-c("mu","pre")
  fit_cont <- rjags::coda.samples(fit_cont, variable.names = params_CONT, n.iter = n.iter, thin = 1)
  out<-do.call(rbind.data.frame, fit_cont)
  # Posterior samples
  mu_CONT<-out$mu
  
  # Difference of posterior distributions by strata
  theta<-mu_TREAT - mu_CONT
  
  # Posterior mean with credible interval
  m<-mean(theta)
  CrI<-as.numeric(quantile(theta, c(za/2,1-za/2)))
  LCrI<-CrI[2] - CrI[1]
  w<-c(m, CrI) 
  

  ######################################################################################################
  ######################################################################################################
  # Trace plot of mean difference
  df<-data.frame(iteration=1:n.iter, theta=theta)
  fp4<-ggplot2::ggplot(df, ggplot2::aes(x=iteration, y=theta)) +
    ggplot2::geom_line(color = "black", linewidth = 0.5, alpha = 0.9, linetype = 1) +
    ggplot2::xlab('Iteration') + 
    ggplot2::ylab('theta') +
    ggplot2::ggtitle('Trace Plot')
  
  #X11(width = 20, height=8)
  #print(fp4)
  ######################################################################################################
  ######################################################################################################
  # Density of posterior samples
  df<-data.frame(x=theta)
  fp5<-ggplot2::ggplot(df, ggplot2::aes(x) ) +
    geom_density(alpha = 0.7, bw=0.05) +
    ggplot2::labs(title = 'Density of posterior samples ', 
         x='theta', y='Density')
  #X11()
  #print(fp5)
  ######################################################################################################
  ######################################################################################################
  # Autocorrelation function for trace plot
  acf_dat <- acf(theta, plot = FALSE)
  df <- with(acf_dat, data.frame(lag, acf))
  ciline <- qnorm(za/2)/sqrt(length(theta))
  
  fp6<-ggplot2::ggplot(data=df, ggplot2::aes(x=lag, y=acf)) + 
    geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_errorbar(ggplot2::aes(x=lag, ymax=acf, ymin=0), width=0) + 
    geom_hline(ggplot2::aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
    geom_hline(ggplot2::aes(yintercept = -ciline), linetype = 2, color = 'darkblue')
  #####################################################################################################
  ######################################################################################################
  # Trace, density, and ACF 
  fp7 <- ggpubr::ggarrange(fp4, fp5, fp6, nrow = 3, ncol = 1)
  
  dev.new()
  print(fp7)
  
  ######################################################################################################
  ######################################################################################################
  # Plot - Incidence proportion for ICA and ECA
  ######################################################################################################
  ######################################################################################################
  label <- c('ICA', 'ECA')
  mean  <- c(R_ICA[1], R_ECA[1]) 
  lower <- c(R_ICA[2], R_ECA[2]) 
  upper <- c(R_ICA[3], R_ECA[3]) 
  
  df <- data.frame(label, mean, lower, upper)
  
  # Reverse the factor level ordering for labels after ggplot2::coord_flip()
  df$label <- factor(df$label, levels=rev(df$label))
  
  fp9 <- ggplot2::ggplot(data=df, ggplot2::aes(x=label, y=mean, ymin=lower, ymax=upper)) +
    ggplot2::geom_pointrange(linewidth = 1) + 
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=0 after flip
    ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
    ggplot2::xlab('') + 
    ggplot2::ylab('Posterior Mean (95% CrI)') +
    ggplot2::ggtitle('Overall event rate for ICA and ECA') +
    ggplot2::theme_bw()  # use a white background
  
  dev.new()
  print(fp9)
  
  ######################################################################################################
  ######################################################################################################
  # Function output
  output<-list(n_1 = n_1, n_0 = n_0, n_11 = n_11, n_10 = n_10, 
               r_ita=R_ITA, r_ica=R_ICA, r_eca=R_ECA,
               A_eff=sum(A*gPPP),
               PPP=PPP, gPPP=gPPP, gamma=gamma, 
               LCrI=LCrI,
               Estimate = w)
  return(output)
  # function ends
}


# ######################################################################################################
# ######################################################################################################
# # Example: Simulate Data
# ######################################################################################################
# ######################################################################################################
# # Parameters
# N0<-3000 # number of external subjects
# N1<-600 # number of current subjects
# 
# p<-10 # number of variables
# rho<-0.1 # correlation between the p variables 
# mu0<-rep(1.1, p) # mean of variables from external data source 
# mu1<-rep(1, p) # mean of variables from current study
# var0<-0.12  # variance of variables from external data source
# var1<-0.1 # variance of variables from current study
# 
# # Covariance matrices
# cov0<-rep(sqrt(var0*rho),p)%*%t(rep(sqrt(var0*rho),p))+diag(p)*var0*(1-rho)
# cov1<-rep(sqrt(var1*rho),p)%*%t(rep(sqrt(var1*rho),p))+diag(p)*var1*(1-rho)
# 
# # Set pseudo-random generator seed
# #set.seed(2*pi)
# 
# # Generate covariates
# x0<-rmvnorm(N0,mu0,cov0)
# x1<-rmvnorm(N1,mu1,cov1)
# 
# # Convert first four variables into binary variables
# x0[,1:4]<-(x0[,1:4]>1)*1
# x1[,1:4]<-(x1[,1:4]>1)*1
# 
# # True coefficients 
# beta<-c(0.2,0.4,0.1,0.3,1,1,1,1.1,1,1)
# 
# # True treatment effect
# delta_treat1=-7
# # True mean difference between internal and external controls
# delta_treat0=0.5
# 
# TRT0<-rep(0, N0)
# TRT1<-rbinom(N1, size=1, prob=2/3)
# 
# # Outcome
# y0 <- 100 + TRT0 + delta_treat0 + x0%*%beta  + rnorm(N0, mean = 0, sd=10)
# y1 <- 100 + TRT1*delta_treat1 + x1%*%beta  + rnorm(N1, mean = 0, sd=5)
# 
# # Create data object including covariates (X), data source indicator (Z), and outcome (Y)
# Y<-c(y0,y1)
# TRT<-c(TRT0, TRT1)
# Z<-c(rep(0,N0),rep(1,N1)) # data source indicator 0: external, 1: current
# X<-rbind(x0,x1)
# 
# # Define data object
# data<-data.frame(X, Z, TRT, Y)
# 
# #####################################################################################################
# ######################################################################################################
# # Example with simulated data
# ######################################################################################################
# ######################################################################################################
# # Variable names for the adjustment
# var.names<-c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10') 
# 
# set.seed(pi)  
# out<-BDB_Outcome_Adj_continuous(Y=data$Y, 
#                                 Z=data$Z, 
#                                 TRT=data$TRT, 
#                                 X = data[,var.names], 
#                                 A = 200, 
#                                 a=1,
#                                 za=0.05,
#                                 n.iter=100)
# 
