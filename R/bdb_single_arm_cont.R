#' Bayesian Dynamic Borrowing for Single-Arm Trial with Continuous Outcome
#'
#' Applies the Bayesian Dynamic Borrowing double-adjustment method for single-arm
#' trials with continuous outcomes. Requires JAGS. The single-arm trial is augmented with external data.
#'
#' @param Y Numeric vector of continuous outcomes.
#' @param Z Numeric vector indicating data source: 1 for current trial data, 0 for external data.
#' @param X Matrix or data frame of baseline covariates (confounders).
#' @param A Numeric. Number of subjects intended to be borrowed from external data source.
#' @param s Integer. Number of strata for propensity score stratification (default 3).
#' @param a Numeric. Tuning parameter for the arctangent elastic function (default 1).
#' @param za Numeric. Significance level for credible intervals (default 0.05 for 95\% CrI).
#' @param n_iter Integer. Number of MCMC iterations (default 5000).
#' @param double_adj Logical. If TRUE (default), both baseline and outcome adjustment; if FALSE, baseline only.
#' @param alpha Numeric. Shape hyperparameter for Gamma prior on precision (default 0.01).
#' @param beta Numeric. Rate hyperparameter for Gamma prior on precision (default 0.01).
#' @param tau Numeric. Precision for Normal prior on mean (default 1e-5).
#' @param show_plots Logical. If TRUE, returns diagnostic plots in \code{plots} element (default FALSE).
#'
#' @return A list with \code{N0_t}, \code{n_s0}, \code{n_s1}, \code{r_s0}, \code{r_s1}, \code{r_sat}, \code{r_eca},
#'   \code{A_eff}, \code{v}, \code{r}, \code{gamma}, \code{PPP}, \code{gPPP}, \code{mean.theta}, \code{var.theta},
#'   \code{CrI}, \code{Width}, \code{theta_s}, \code{theta_o}, \code{theta} (same as \code{theta_o}),
#'   and optionally \code{plots} when \code{show_plots = TRUE}.
#'
#' @references
#' Farjat, A., Ji, Y., Kaiser, A., Lu, C., Pap, A., Potts, J., Wang, M. A Comprehensive Bayesian Double-Adjustment
#' Approach to Dynamic Borrowing of External Data. Manuscript in progress.
#' Wang, C., Li, H., Chen, W.C., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2019. Journal of biopharmaceutical statistics, 29(5), pp.731-748.
#'
#' @seealso \code{\link{bdb_single_arm_bin}}, \code{\link{bdb_two_arms_bin}}, \code{\link{bdb_two_arms_cont}}, \code{\link{elastic_function}}
#'
#' @examples
#' \dontrun{
#' # Requires rjags and JAGS.
#' # See vignette("Introduction to BayesDDB", package = "BayesDDB") for full examples.
#' }
#' @export
bdb_single_arm_cont <- function(Y, Z, X, A, s = 3, a = 1, za = 0.05, n_iter = 5000,
                                double_adj = TRUE, alpha = 0.01, beta = 0.01, tau = 1e-5, show_plots = FALSE) {
  if (!requireNamespace("rjags", quietly = TRUE))
    stop("Package 'rjags' is required. Install it and ensure JAGS is installed.", call. = FALSE)
  for (pkg in c("dplyr", "overlapping", "ggplot2", "ggpubr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Please install it.", call. = FALSE)
  }
  doubleadj <- double_adj
  plots <- show_plots
  n.iter <- n_iter
  model_string_PPP <- "
  model
  {
    for (i in 1:n){
      y[i] ~ dnorm(m,pre)
   }
  m ~ dnorm(0,tau)
  pre ~ dgamma(alpha,beta)
  new ~ dnorm(m,pre*N1)
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
    mu ~ dnorm(0, tau)
    pre ~ dgamma(alpha, beta)
  }
  ", file = modelfile)

  # Define data object
  data<-data.frame(Y, Z, X)
  
  ######################################################################################################
  # Step 1: Modeling 
  ######################################################################################################
  # Define exposure mechanism
  mylogit<-glm(Z ~ .-Y, data = data, family = stats::binomial())
  
  ######################################################################################################
  # Step 2: Propensity Score Calculation 
  ######################################################################################################
  prd<-predict(mylogit, type='response')
  data<-data.frame(data, prd)
  prd1<-prd[Z == 1] # propensity scores of current patients
  
  ######################################################################################################
  # Step 3: Trimming 
  ######################################################################################################
  # Define data1 and restrict propensity scores to the range of values of the 
  # current patients, this is just trimming
  data1<-data[dplyr::between(prd, min(prd1), max(prd1)),]
  # Propensity score of patients in the trial
  prd1<-data1$prd[data1$Z==1]
  # Propensity scores of external patients
  prd0<-data1$prd[data1$Z==0]
  # Overall number of subjects after trimming
  n<-nrow(data1) 
  # Number of current subjects
  N1<-sum(Z)
  # Number of external subjects after trimming
  n0<-n-N1
  
  ######################################################################################################
  # Step 4: Stratification 
  ######################################################################################################
  # s, number of strata
  if (s==1) {q=c(0,1)
  } else { 
    q<-quantile(prd1, seq(1:(s-1))/s)
    q<-c(0, q, 1) 
  }
  # Define matrix strat to identify indices of each stratum 
  strat<-matrix(data=rep(FALSE, length(data1$prd)*s), nrow = s, ncol = length(data1$prd))
  for (i in 1:s){ strat[i,]<-dplyr::between(data1$prd, q[i], q[i+1]) }
  
  # Define variable s_ind
  data1$stratum<-NA
  for(i in 1:s) {data1$stratum[strat[i,]==1]<-i}
  
  ######################################################################################################
  # Step 5: Overlapping 
  ######################################################################################################
  
  v<-rep(0,s)
  for (i in 1:s){
    # overlapping probability of stratum i
    temp_lst<-list(P0_s=data1$prd[data1$Z==0 & data1$stratum==i], 
                   P1_s=data1$prd[data1$Z==1 & data1$stratum==i])
    v[i]<-as.numeric(overlapping::ovmult(x=temp_lst)) 
  }
  
  ######################################################################################################
  # Step 6; Weighting 
  ######################################################################################################
  # Initialize vector of normalized overlapping coefficients 
  r<-rep(0,s)
  for (i in 1:s){ r[i]<-v[i]/sum(v) }
  
  ######################################################################################################
  # Step 7: Calibration by stratum
  ######################################################################################################
  # Estimate predictive posterior probability by stratum
  e.m<-rep(0,s) # observed mean in SAT by stratum
  for(i in 1:s){ e.m[i]<-mean(data1$Y[data1$Z==1 & data1$stratum==i]) }
  
  PPP<-rep(0,s)
  for(i in 1:s){
    # jags model  
    model_PPP<-textConnection(model_string_PPP)
    # Initial values of parameters
    inits<-list(m=e.m[i])
    # Data list as input to JAGS
    n0_s<-sum(data1$Z==0 & data1$stratum==i)
    y0_s<-data1$Y[data1$Z==0 & data1$stratum==i]
    N1_s<-sum(data1$Z==1 & data1$stratum==i)
    #if(N1_s==0) {N1_s=1} # otherwise the model does not run
    dataList<-list("n"=n0_s,"y"=y0_s,"N1"=N1_s, "alpha"=alpha, "beta"=beta, "tau"=tau) 
    
    # Run MCMC by JAGS (sampling from external data assuming precision from current data)
    jagsfit <- rjags::jags.model(model_PPP, data=dataList, inits=inits, n.chains=1, n.adapt=10000) 
    # Draw MCMC samples
    params<-c("m","new") # Parameters to be monitored in MCMC
    out <- rjags::coda.samples(jagsfit, params, n.iter=n.iter)
    out<-do.call(rbind.data.frame, out)
    # Calculate posterior predictive probability
    PPP[i]<-mean(out$new<e.m[i])
  }
  
  ######################################################################################################
  # Step 8: Mapping 
  ######################################################################################################
  # a: value of the tuning parameter in the arc-tangent elastic function
  gPPP<-atan(a*sin(PPP*pi))/atan(a) # further discounting by the elastic function 
  
  # In case of baseline adjustment only
  if(doubleadj==FALSE) {gPPP=rep(1,s)} 
  
  ######################################################################################################
  # Step 9: Discounting 
  ######################################################################################################
  n_s0<-rep(0,s)
  for (i in 1:s){ n_s0[i]<-sum(data1$Z==0 & data1$stratum==i) }
  
  n_s1<-rep(0,s)
  for (i in 1:s){ n_s1[i]<-sum(data1$Z==1 & data1$stratum==i) }
  
  # Discount parameter
  gamma<-rep(0,s)
  for (i in 1:s){ gamma[i]<-min(1,A*r[i]/n_s0[i]) }
  
  # Mean and 95% CI for SAT
  R_SAT<-mean_CI(data1$Y[data1$Z==1], za)
  # Mean and 95% CI for ECA
  R_ECA<-mean_CI(data1$Y[data1$Z==0], za)
  # Mean and 95% CI for ECA by stratum
  r_s0<-matrix(data=NA, nrow=s, ncol=3)
  for(i in 1:s){r_s0[i,]<-mean_CI(data1$Y[data1$Z==0 & data1$stratum==i], za)}
  # Mean and 95% CI for SAT by stratum
  r_s1<-matrix(data=NA, nrow=s, ncol=3)
  for(i in 1:s){r_s1[i,]<-mean_CI(data1$Y[data1$Z==1 & data1$stratum==i], za)}
  
  ######################################################################################################
  # Steps 10-11 Analysis/Summary 
  ######################################################################################################
  # Initialize post.mean.samples object 
  #post.mean.samples<-matrix(data=rep(0, s*n.iter), nrow=s, ncol=n.iter)
  # Initialize theta and theta_hat_s objects
  theta<-matrix(rep(0, s*n.iter), nrow=s, ncol=n.iter)
  theta_hat_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  # Obtain posterior mean samples by strata
  for (i in 1:s){ 
    # Internal and external controls
    # Data list
    n_<-sum(data1$stratum==i)
    y_<-data1$Y[data1$stratum==i]
    gamma_<-data1$Z[data1$stratum==i] + (data1$Z[data1$stratum==i]==0)*gamma[i]*gPPP[i]
    zeros_<-rep(0,n_) 
    data_CONT<-list(n=n_, y=y_, gamma=gamma_, zeros=zeros_, "alpha"=alpha, "beta"=beta, "tau"=tau)
    # Initial values
    mu_<-mean(y_)
    pre_<-1/var(y_)
    inits_CONT<-function(){ list( mu=mu_, pre=pre_) }
    # Fit model 
    fit_cont_s <- rjags::jags.model(file = modelfile, data = data_CONT, inits = inits_CONT, n.chains = 1, n.adapt = 10000)
    update(fit_cont_s, n.iter)
    # Parameters to monitor during MCMC
    params_CONT<-c("mu","pre")
    fit_cont_s <- rjags::coda.samples(fit_cont_s, variable.names = params_CONT, n.iter = n.iter, thin = 1)
    out<-do.call(rbind.data.frame, fit_cont_s)
    # Posterior samples
    theta[i,]<-out$mu
    # Posterior mean by stratum and 95% CrI
    theta_hat_s[,i]<-c(mean(theta[i,]), as.numeric(quantile(theta[i,],c(za/2,1-za/2)))) 
  }
  
  post.dist<-colMeans(theta)
  
  post.mean<-mean(post.dist) # posterior mean
  
  post.var<-var(post.dist)
  
  CrI<-as.numeric(quantile(post.dist,c(za/2,1-za/2))) # credible interval
  
  Width<-CrI[2]-CrI[1] # width of credible interval
  
  # Overall estimate and CrI
  theta_o<-c(post.mean, CrI)
  
  plots_list <- NULL
  if (plots == TRUE) {
    temp <- 'Stratum 1'
    if (s > 1) { for (i in 2:s) { temp <- c(temp, paste('Stratum', i)) } }
    temp <- c(temp, 'Overall')
    label <- temp
    mean  <- c(theta_hat_s[1,], theta_o[1])
    lower <- c(theta_hat_s[2,], theta_o[2])
    upper <- c(theta_hat_s[3,], theta_o[3])
    df <- data.frame(label, mean, lower, upper)
    df$label <- factor(df$label, levels = rev(df$label))
    fp1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = label, y = mean, ymin = lower, ymax = upper)) +
      ggplot2::geom_pointrange(linewidth = 1) + ggplot2::coord_flip() +
      ggplot2::geom_text(x = 1.2, y = theta_o[1], label = paste(round(theta_o[1], 2), ' [', round(theta_o[2], 2), ', ', round(theta_o[3], 2), ']', sep = '')) +
      ggplot2::xlab('') + ggplot2::ylab('Posterior Mean (95% CrI)') + ggplot2::ggtitle('Mean estimate by Stratum and Overall') + ggplot2::theme_bw()
    df <- data.frame(x = data1$prd, Control = c(rep('Current', sum(data1$Z)), rep('External', sum(data1$Z == 0))))
    fp2 <- ggplot2::ggplot(df, ggplot2::aes(x = x, fill = Control)) + ggplot2::geom_density(alpha = 0.7, bw = 0.05) +
      ggplot2::labs(title = 'Propensity Score Distribution by Group', x = 'Propensity Score', y = 'Density', fill = 'Group')
    df <- data.frame(prd = data1$prd, stratum = data1$stratum, Control = c(rep('Current', sum(data1$Z)), rep('External', sum(data1$Z == 0))))
    fp3 <- ggplot2::ggplot(df, ggplot2::aes(x = prd, fill = Control)) + ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap(~ stratum, scales = 'free') +
      ggplot2::labs(title = "Propensity Score Distribution by Group and Stratum", x = "Propensity Score", y = "Density", fill = "Group") + ggplot2::theme_minimal()
    df <- data.frame(iteration = 1:n.iter, theta = post.dist)
    fp4 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = theta)) + ggplot2::geom_line(color = "black", linewidth = 0.5, alpha = 0.9, linetype = 1) + ggplot2::xlab('Iteration') + ggplot2::ylab('theta') + ggplot2::ggtitle('Trace Plot')
    df <- data.frame(x = post.dist)
    fp5 <- ggplot2::ggplot(df, ggplot2::aes(x)) + ggplot2::geom_density(alpha = 0.7, bw = 0.05) + ggplot2::labs(title = 'Density of posterior samples', x = 'theta', y = 'Density')
    acf_dat <- acf(post.dist, plot = FALSE)
    df <- with(acf_dat, data.frame(lag, acf))
    ciline <- qnorm(za/2) / sqrt(length(post.dist))
    fp6 <- ggplot2::ggplot(data = df, ggplot2::aes(x = lag, y = acf)) + geom_hline(ggplot2::aes(yintercept = 0)) + ggplot2::geom_errorbar(ggplot2::aes(x = lag, ymax = acf, ymin = 0), width = 0) + geom_hline(ggplot2::aes(yintercept = ciline), linetype = 2, color = 'darkblue') + geom_hline(ggplot2::aes(yintercept = -ciline), linetype = 2, color = 'darkblue')
    fp7 <- ggpubr::ggarrange(fp4, fp5, fp6, nrow = 3, ncol = 1)
    label <- c('ECA', 'SAT')
    mean  <- c(R_ECA[1], R_SAT[1])
    lower <- c(R_ECA[2], R_SAT[2])
    upper <- c(R_ECA[3], R_SAT[3])
    group <- c('ECA', 'SAT')
    df <- data.frame(label, mean, lower, upper, group)
    df$label <- factor(df$label, levels = rev(df$label))
    df$group <- factor(df$group, levels = rev(df$group))
    fp9 <- ggplot2::ggplot(data = df, ggplot2::aes(x = label, y = mean, ymin = lower, ymax = upper, fill = group)) +
      ggplot2::geom_linerange(linewidth = 5, position = ggplot2::position_dodge(width = 0.25), colour = "lightgrey") +
      ggplot2::geom_point(size = 3, shape = 21, stroke = 0.25, position = ggplot2::position_dodge(width = 0.25)) +
      ggplot2::coord_flip() + ggplot2::xlab('') + ggplot2::ylab('Posterior Mean (95% CrI)') + ggplot2::ggtitle('Overall sample mean for ECA and SAT') + ggplot2::theme_bw()
    temp <- 'Stratum 1'
    if (s > 1) { for (i in 2:s) { temp <- c(temp, paste('Stratum', i)) } }
    temp <- c(temp, temp)
    label <- temp
    group <- c(rep('ECA', s), rep('SAT', s))
    rate  <- c(r_s0[,1], r_s1[,1])
    lower <- c(r_s0[,2], r_s1[,2])
    upper <- c(r_s0[,3], r_s1[,3])
    df <- data.frame(group, label, rate, lower, upper)
    df$label <- factor(df$label, levels = rev(df$label)[1:s])
    df$group <- factor(df$group, levels = c('SAT', 'ECA'))
    fp10 <- ggplot2::ggplot(data = df, ggplot2::aes(x = label, y = rate, ymin = lower, ymax = upper, fill = group)) +
      ggplot2::geom_linerange(linewidth = 5, position = ggplot2::position_dodge(width = 0.25), colour = "lightgrey") +
      ggplot2::geom_point(size = 3, shape = 21, stroke = 0.25, position = ggplot2::position_dodge(width = 0.25)) +
      ggplot2::scale_x_discrete(name = " ") + ggplot2::scale_y_continuous(name = "Sample mean [95% CI]") + ggplot2::coord_flip() +
      ggplot2::ggtitle('Sample mean for SAT and ECA by stratum') + ggplot2::theme_bw()
    plots_list <- list(forest = fp1, ps_group = fp2, ps_stratum = fp3, mcmc = fp7, overall_rates = fp9, rates_by_stratum = fp10)
  }
  output <- list(N0_t = n0, n_s0 = n_s0, n_s1 = n_s1, r_s0 = r_s0, r_s1 = r_s1,
                 r_sat = R_SAT, r_eca = R_ECA,
                 A_eff = A * sum(r * gPPP), v = v, r = r, gamma = gamma, PPP = PPP, gPPP = gPPP,
                 mean.theta = post.mean, var.theta = post.var, CrI = CrI, Width = Width,
                 theta_s = theta_hat_s, theta_o = theta_o, theta = theta_o)
  if (!is.null(plots_list)) output$plots <- plots_list
  output
}

