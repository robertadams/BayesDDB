#' Outcome-adjustment-only helper for binary two-arm case (internal)
#'
#' Applies outcome adjustment for Bayesian Dynamic Borrowing in the two-arm binary setting.
#' Used internally; not called by the main \code{bdb_*} API.
#'
#' @param Y Numeric vector of binary outcomes.
#' @param Z Numeric vector: 1 = current trial, 0 = external.
#' @param TRT Numeric vector: 1 = treatment, 0 = control.
#' @param X Matrix or data frame of baseline covariates.
#' @param A Numeric. Intended number of subjects to borrow.
#' @param a Numeric. Tuning parameter for elastic function.
#' @param s Integer. Number of strata.
#' @param alpha Numeric. Beta prior shape1.
#' @param beta Numeric. Beta prior shape2.
#' @param za Numeric. Significance level for credible intervals.
#' @param K Integer. Number of posterior samples.
#' @return List with estimate, credible interval, and diagnostic quantities.
#' @keywords internal
#' @noRd
BDB_Outcome_Adj_binary <- function(Y, Z, TRT, X, A, a, s, alpha, beta, za, K) {
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
  # alpha: shape1 hyper-parameter from prior beta distribution of binomial distribution success parameter p 
  # beta : shape2 hyper-parameter from prior beta distribution of binomial distribution success parameter p
  # za: confidence level value
  # K: number of samples from posterior distribution
  #
  # OUTPUT
  # n_1: number of current subjects  
  # n_0: number of external subjects
  # n_11: current subjects in treatment arm 
  # n_10: number of internal control subjects 
  # ne_11: number of events in current treatment arm 
  # ne_10: number of events in current control arm 
  # ne_0: number of events in external data source
  # R_ITA: overall event rate and 95% CrI of internal treatment arm
  # R_ICA: overall event rate and 95% CrI of internal control arm
  # R_ECA: overall event rate and 95% CrI of external control arm
  # A_eff: Overall effective number of subjects borrowed from external sources
  # PPP: posterior predictive probability
  # gPPP: transformed posterior predictive probability used for further discounting
  # gamma: s-by-1 vector of power discounting parameter
  # LCrI: length of credible interval from Estimate
  # Estimate: 1-by-3 vector of treatment effect and corresponding 95% CrI. 
  # The treatment effect is defined as the difference of event rates between treatment 
  # and control arms 
  #
  # Plots
  # Figure 1: Trace plot, density and ACF from posterior samples of of rate difference
  # Figure 2: Overall event rate and 95% CrI for internal (ICA) and external control arm (ECA)
  ######################################################################################################
  # Function begins
  # Define data object
  data<-data.frame(Y, Z, TRT, X)

  ######################################################################################################
  ######################################################################################################
  # Step 7: Calibration 
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
  # Number of events for current treatment 
  ne_11<-sum(data$Y[data$Z==1 & data$TRT==1])
  # Number of events for internal control
  ne_10<-sum(data$Y[data$Z==1 & data$TRT==0]) 
  # Number of events for external data
  ne_0<- sum(data$Y[data$Z==0]) 
  
  # Estimate predictive posterior probability
  e.m<-ne_10/n_10 # observed mean of current control
  
  # Number of random samples
  N<-10000
  # beta-binomial posterior for external data 
  prb <- rbeta(N, shape1 = alpha + ne_0, shape2 = beta + n_0 - ne_0) 
  
  # Posterior predictive probability (sampling from current study with probabilities from external data and 
  # compare against current study)
  PPP<-mean(rbinom(N, n_10, prb)/n_10 > e.m)
  
  ######################################################################################################
  ######################################################################################################
  # Step 8: Mapping 
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
  # Sampling from beta-binomial distribution with parameters shape parameters alpha and beta
  thetaT<-rbeta(K, alpha + ne_11, beta + n_11-ne_11)
  # Current and external controls
  # Sampling from beta-binomial distribution with parameters shape parameters alpha and beta
  thetaC<-rbeta(K, alpha + ne_10 + ne_0*gamma*gPPP, beta + (n_10-ne_10) + (n_0-ne_0)*gamma*gPPP)
  # Treatment effect (Treatment - Controls)
  theta<-thetaT - thetaC
  
  # Posterior mean with credible interval
  m<-mean(theta)
  CrI<-as.numeric(quantile(theta, c(za/2,1-za/2)))
  LCrI<-CrI[2] - CrI[1]
  w<-c(m, CrI) 

  
  ######################################################################################################
  ######################################################################################################
  # Current treatment event rate 
  ######################################################################################################
  ######################################################################################################
  # Event rate and CI (Clopper-Person)
  R_ITA<-c(ne_11/n_11 , qbeta(za/2, ne_11, n_11 - ne_11 +1) , qbeta(1-za/2, ne_11 + 1, n_11 - ne_11))
  
  ######################################################################################################
  ######################################################################################################
  # ICA event rate 
  ######################################################################################################
  ######################################################################################################
  # Overall event rate from ICA and confidence interval from beta distribution (Clopper-Person interval)
  R_ICA<-c(ne_10/n_10 , qbeta(za/2, ne_10, n_10 - ne_10 + 1) , qbeta(1-za/2, ne_10 + 1, n_10 - ne_10))
  
  ######################################################################################################
  ######################################################################################################
  # ECA event rate 
  ######################################################################################################
  ######################################################################################################
  # Overall event rate from ECA and confidence interval from beta distribution (Clopper-Person interval)
  R_ECA<-c(ne_0/n_0, qbeta(za/2, ne_0, n_0 - ne_0 +1) , qbeta(1-za/2, ne_0 + 1, n_0 - ne_0))
  
  ######################################################################################################
  ######################################################################################################
  # Trace plot of mean difference
  df<-data.frame(iteration=1:K, theta=theta)
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
               ne_11=ne_11, ne_10=ne_10, ne_0=ne_0,
               r_ita=R_ITA, r_ica=R_ICA, r_eca=R_ECA,
               A_eff=sum(A*gPPP),
               PPP=PPP, gPPP=gPPP, gamma=gamma, 
               LCrI=LCrI,
               Estimate = w)
  return(output)
}


