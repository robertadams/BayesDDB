#' Bayesian Dynamic Borrowing for Single-Arm Trial with Binary Outcome
#'
#' Applies the Bayesian Dynamic Borrowing double-adjustment method for single-arm
#' trials with binary outcomes. The single-arm trial is augmented with external
#' data. No treatment variable is used.
#'
#' @param Y Numeric vector of binary outcomes (0s and 1s).
#' @param Z Numeric vector indicating data source: 1 for current trial data,
#'   0 for external data.
#' @param X Matrix or data frame of baseline covariates (confounders) with
#'   n rows and p columns.
#' @param A Numeric. Number of subjects intended to be borrowed from external
#'   data source.
#' @param s Integer. Number of strata for propensity score stratification (default 3).
#' @param a Numeric. Tuning parameter for the arctangent elastic function (must be >= 1, default 1).
#' @param alpha Numeric. Shape1 hyperparameter for the Beta prior (default 0.5).
#' @param beta Numeric. Shape2 hyperparameter for the Beta prior (default 0.5).
#' @param za Numeric. Significance level for credible intervals (default 0.05 for 95\% CrI).
#' @param K Integer. Number of samples from the posterior distribution (default 10000).
#' @param double_adj Logical. If TRUE (default), both baseline and outcome adjustment. If FALSE, baseline only.
#' @param show_plots Logical. If TRUE, returns diagnostic plots in \code{plots} element (default FALSE).
#'
#' @return A list containing \code{n_s1}, \code{n_s0}, \code{ne_s1}, \code{ne_s0},
#'   \code{r_s1}, \code{r_s0}, \code{r_sat}, \code{r_eca}, \code{v}, \code{r},
#'   \code{A_eff}, \code{PPP}, \code{gPPP}, \code{gamma}, \code{mean.theta},
#'   \code{var.theta}, \code{CrI}, \code{Width}, \code{theta_s}, \code{theta_o},
#'   and optionally \code{plots} (list of ggplot objects) when \code{show_plots = TRUE}.
#'
#' @references
#' Farjat, A., Ji, Y., Kaiser, A., Lu, C., Pap, A., Potts, J., Wang, M.
#' A Comprehensive Bayesian Double-Adjustment Approach to Dynamic Borrowing of
#' External Data. Manuscript in progress.
#'
#' Wang, C., Li, H., Chen, W.C., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2019.
#' Propensity score-integrated power prior approach for incorporating real-world
#' evidence in single-arm clinical studies.
#' Journal of biopharmaceutical statistics, 29(5), pp.731-748.
#'
#' @seealso \code{\link{bdb_single_arm_cont}}, \code{\link{bdb_two_arms_bin}},
#'   \code{\link{bdb_two_arms_cont}}, \code{\link{elastic_function}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n0 <- 50
#' n1 <- 30
#' X <- matrix(rnorm((n0 + n1) * 2), ncol = 2)
#' Z <- c(rep(0, n0), rep(1, n1))
#' Y <- c(rbinom(n0, 1, 0.4), rbinom(n1, 1, 0.35))
#' bdb_single_arm_bin(Y = Y, Z = Z, X = X, A = 10, s = 2, a = 1, K = 500, show_plots = FALSE)
#' }
#' @export
bdb_single_arm_bin <- function(Y, Z, X, A, s = 3, a = 1, alpha = 0.5, beta = 0.5,
                               za = 0.05, K = 10000, double_adj = TRUE, show_plots = FALSE) {
  for (pkg in c("dplyr", "overlapping", "ggplot2", "ggpubr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. Please install it.", call. = FALSE)
  }
  doubleadj <- double_adj
  plots <- show_plots
  data <- data.frame(Y, Z, X)
  ######################################################################################################
  # Step 1: Model exposure mechanism
  ######################################################################################################
  # Define exposure mechanism
  mylogit <- glm(Z ~ . -Y, family = stats::binomial(), data = data)
  ######################################################################################################
  ######################################################################################################
  # Step 2: Propensity Score Calculation 
  ######################################################################################################
  ######################################################################################################
  # Propensity scores
  prd<-predict(mylogit, type='response')
  # Add propensity scores to data object
  data<-data.frame(data, prd)
  # Propensity scores of current patients (patients in the trial, RCT: Z=1, ECA: Z=0)
  prd1<-prd[Z == 1] 
  
  ######################################################################################################
  ######################################################################################################
  # Step 3: Trimming 
  ######################################################################################################
  ######################################################################################################
  # Define data1 and restrict propensity scores to the range of values of the 
  # current patients, this is just trimming
  data1<-data[dplyr::between(prd, min(prd1), max(prd1)),]
  # Sort data by Z
  data1<-data1[order(data1$Z, decreasing = TRUE),] 
  # Trimmed Propensity scores
  # Propensity score of patients in the trial
  prd1<-data1$prd[data1$Z==1]
  # Propensity scores of external patients
  prd0<-data1$prd[data1$Z==0]
  
  ######################################################################################################
  ######################################################################################################
  # Step 4: Stratification 
  ######################################################################################################
  ######################################################################################################
  # s,  number of strata
  if (s==1) {q=c(0,1)
  } else { 
    q<-quantile(prd1, seq(1:(s-1))/s)
    q<-c(0, q, 1) 
  }
 
  # define matrix strat to identify indices of each stratum 
  strat<-matrix(data=rep(FALSE, length(data1$prd)*s), nrow = s, ncol = length(data1$prd))
  for (i in 1:s){
    strat[i,]<-dplyr::between(data1$prd, q[i], q[i+1])
  }
  
  # Define variable stratum 
  data1$stratum<-NA
  for(i in 1:s) {data1$stratum[strat[i,]==1]<-i}
  
  # Current trial
  # n_s1, number subjects in the current trial by stratum
  n_s1<-rep(0,s)
  for (i in 1:s){ n_s1[i]<-sum(data1$Z==1 & data1$stratum==i) }
  # n_1, overall number of subjects in the current trial
  n_1<-sum(n_s1)
  # ne_s1, number of events in current trial by stratum 
  ne_s1<-rep(0,s)
  for(i in 1:s){ ne_s1[i]<-sum(data1$Y[data1$Z==1 & data1$stratum==i]) }
  # ne_1, overall number of events in the current trial 
  ne_1<-sum(ne_s1)
  
  # External Control
  # n_s0, number of subject in external control by stratum
  n_s0<-rep(0,s)
  for (i in 1:s){ n_s0[i]<-sum(data1$Z==0 & data1$stratum==i) }
  # n_0, overall number of subjects in the external control
  n_0<-sum(n_s0)
  # ne_s0, number of event for external data by stratum
  ne_s0<-rep(0,s)
  for(i in 1:s){ ne_s0[i]<- sum(data1$Y[data1$Z==0 & data1$stratum==i]) }
  # ne_0, overall number of events in the external in control 
  ne_0<-sum(ne_s0)
  
  ######################################################################################################
  ######################################################################################################
  # Step 5: Overlapping coefficient 
  ######################################################################################################
  ######################################################################################################
  # Calculate vector of overlapping coefficients
  v<-rep(0,s)
  for (i in 1:s){
    # overlapping probability of stratum i
    temp_lst<-list(P0_s=data1$prd[data1$Z==0 & data1$stratum==i], 
                   P1_s=data1$prd[data1$Z==1 & data1$stratum==i])
    v[i]<-as.numeric(overlapping::ovmult(x=temp_lst)) 
  }
  
  ######################################################################################################
  ######################################################################################################
  # Step 6: Weighting 
  ######################################################################################################
  ######################################################################################################
  # Initiate vector of normalized overlapping coefficients 
  r<-rep(0,s)
  for (i in 1:s){
    r[i]<-v[i]/sum(v)
  }
  
  ######################################################################################################
  ######################################################################################################
  # Step 7: Calibration 
  ######################################################################################################
  ######################################################################################################
  # Observed mean in current trial by stratum
  e.m<-ne_s1/n_s1 
  
  # Number of random samples
  N<-10000
  # beta-binomial posterior for external data by stratum
  prb<-matrix(0, nrow = s, ncol=N)
  for(i in 1:s){ prb[i,] <- rbeta(N, shape1 = alpha + ne_s0[i], shape2 = beta + n_s0[i] - ne_s0[i]) }
  
  # Posterior predictive probability (sampling from current study with probabilities from external data and 
  # compare against current study)
  PPP<-rep(0, s)
  for(i in 1:s){ PPP[i]<-mean(rbinom(N, n_s1[i], prb[i,])/n_s1[i] > e.m[i]) }
  
  ######################################################################################################
  ######################################################################################################
  # Step 8: Mapping 
  ######################################################################################################
  ######################################################################################################
  # a: value of the tuning parameter in the arctangent elastic function
  
  # Further discounting by the elastic function 
  gPPP<-atan(a*sin(PPP*pi))/atan(a) # mapped predictive posterior probability
  
  # In case of baseline adjustment only
  if(doubleadj==FALSE) {gPPP=rep(1,s)} 
  
  ######################################################################################################
  ######################################################################################################
  # Step 9: Discounting 
  ######################################################################################################
  ######################################################################################################
  # A: number of subjects intended to be borrowed from external source
  # Define discounting factor
  gamma<-rep(0,s)
  for (i in 1:s){ gamma[i]<-min(1,A*r[i]/n_s0[i]) }
  
  # Number of subjects actually borrowed from external source by stratum
  A_borrow<-rep(0,s)
  for (i in 1:s){ A_borrow[i]<-A*r[i] }
  
  ######################################################################################################
  ######################################################################################################
  # Steps 10-11: Analysis/Summary 
  ######################################################################################################
  ######################################################################################################
  
  ######################################################################################################
  ######################################################################################################
  # Both baseline and outcome adjustment 
  ######################################################################################################
  ######################################################################################################
  # Initialize matrices for later use
  theta<-matrix(rep(0, s*K), nrow=s, ncol=K) 
  theta_s<-matrix(rep(0, 3*s), nrow=3, ncol=s) 
  for (i in 1:s){
    # Single-arm trial with discounted external data
    theta[i,]<-rbeta(K, alpha + ne_s1[i] + ne_s0[i]*gamma[i]*gPPP[i], beta + (n_s1[i]-ne_s1[i]) + (n_s0[i]-ne_s0[i])*gamma[i]*gPPP[i])
    # Posterior mean by stratum and 95% CrI
    theta_s[,i]<-c(mean(theta[i,]), as.numeric(quantile(theta[i,],c(za/2,1-za/2))))
  }
  # Calculate mean across strata (posterior distribution)
  theta.m<-colMeans(theta)
  post.mean<-mean(theta.m)
  post.var<-var(theta.m)
  # Credible Interval
  CrI<-as.numeric(quantile(theta.m, c(za/2,1-za/2)))
  # Length of credible interval 
  Width<-CrI[2] - CrI[1]
  # Overall posterior mean with credible interval
  theta_o<-c(post.mean, CrI)
  
  ######################################################################################################
  ######################################################################################################
  # Current treatment event rate 
  ######################################################################################################
  ######################################################################################################
  # Event rate by stratum and CI (Clopper-Person)
  r_s1<-t(rbind(ne_s1/n_s1 , qbeta(za/2, ne_s1, n_s1-ne_s1+1) , qbeta(1-za/2, ne_s1+1, n_s1-ne_s1)))
  
  # Overall event rate from SAT and confidence interval from beta distribution (Clopper-Person interval)
  R_SAT<-c(ne_1/n_1, qbeta(za/2, ne_1, n_1 - ne_1 +1) , qbeta(1-za/2, ne_1+1, n_1-ne_1) )
  
  ######################################################################################################
  ######################################################################################################
  # ECA event rate 
  ######################################################################################################
  ######################################################################################################
  # Event rate by stratum and CrI
  r_s0<-t(rbind(ne_s0/n_s0 , qbeta(za/2, ne_s0, n_s0-ne_s0+1) , qbeta(1-za/2, ne_s0+1, n_s0-ne_s0)))
  
  # Overall event rate from ECA and confidence interval from beta distribution (Clopper-Person interval)
  R_ECA<-c(ne_0/n_0, qbeta(za/2, ne_0, n_0 - ne_0 +1) , qbeta(1-za/2, ne_0+1, n_0-ne_0) )
  
  plots_list <- NULL
  if (plots == TRUE) {
    temp <- 'Stratum 1'
    if (s > 1) { for (i in 2:s) { temp <- c(temp, paste('Stratum', i)) } }
    temp <- c(temp, 'Overall')
    label <- temp
    mean  <- c(theta_s[1,], theta_o[1])
    lower <- c(theta_s[2,], theta_o[2])
    upper <- c(theta_s[3,], theta_o[3])
    df <- data.frame(label, mean, lower, upper)
    df$label <- factor(df$label, levels = rev(df$label))
    fp1 <- ggplot2::ggplot(data = df, ggplot2::aes(x = label, y = mean, ymin = lower, ymax = upper)) +
      ggplot2::geom_pointrange(linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::coord_flip() +
      ggplot2::geom_text(x = 1.2, y = theta_o[1], label = paste(round(theta_o[1], 5), ' [', round(theta_o[2], 5), ', ', round(theta_o[3], 5), ']', sep = '')) +
      ggplot2::xlab('') + ggplot2::ylab('Posterior Mean (95% CrI)') + ggplot2::ggtitle('Event rate by Stratum and Overall') + ggplot2::theme_bw()
    df <- data.frame(x = data1$prd, Control = c(rep('Current', sum(data1$Z)), rep('External', sum(data1$Z == 0))))
    fp2 <- ggplot2::ggplot(df, ggplot2::aes(x = x, fill = Control)) +
      ggplot2::geom_density(alpha = 0.7, bw = 0.05) +
      ggplot2::labs(title = 'Propensity Score Distribution by Group', x = 'Propensity Score', y = 'Density', fill = 'Group')
    df <- data.frame(prd = data1$prd, stratum = data1$stratum, Control = c(rep('Current', sum(data1$Z)), rep('External', sum(data1$Z == 0))))
    fp3 <- ggplot2::ggplot(df, ggplot2::aes(x = prd, fill = Control)) +
      ggplot2::geom_density(alpha = 0.5) + ggplot2::facet_wrap(~ stratum, scales = 'free') +
      ggplot2::labs(title = "Propensity Score Distribution by Group and Stratum", x = "Propensity Score", y = "Density", fill = "Group") + ggplot2::theme_minimal()
    df <- data.frame(iteration = 1:K, theta = theta.m)
    fp4 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = theta)) + ggplot2::geom_line(color = "black", linewidth = 0.5, alpha = 0.9, linetype = 1) + ggplot2::xlab('Iteration') + ggplot2::ylab('theta') + ggplot2::ggtitle('Trace Plot')
    df <- data.frame(x = theta.m)
    fp5 <- ggplot2::ggplot(df, ggplot2::aes(x)) + ggplot2::geom_density(alpha = 0.7, bw = 0.05) + ggplot2::labs(title = 'Density of posterior samples', x = 'theta', y = 'Density')
    acf_dat <- acf(theta.m, plot = FALSE)
    df <- with(acf_dat, data.frame(lag, acf))
    ciline <- qnorm(za/2) / sqrt(length(theta.m))
    fp6 <- ggplot2::ggplot(data = df, ggplot2::aes(x = lag, y = acf)) + ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) + ggplot2::geom_errorbar(ggplot2::aes(x = lag, ymax = acf, ymin = 0), width = 0) + ggplot2::geom_hline(ggplot2::aes(yintercept = ciline), linetype = 2, color = 'darkblue') + ggplot2::geom_hline(ggplot2::aes(yintercept = -ciline), linetype = 2, color = 'darkblue')
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
      ggplot2::coord_flip() + ggplot2::xlab('') + ggplot2::ylab('Posterior Mean (95% CrI)') + ggplot2::ggtitle('Overall event rate for ECA and SAT') + ggplot2::theme_bw()
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
      ggplot2::xlab('') + ggplot2::ylab('Event rate [95% CrI]') + ggplot2::coord_flip() + ggplot2::ggtitle('Event rate for ECA and SAT by stratum') + ggplot2::theme_bw()
    plots_list <- list(forest = fp1, ps_group = fp2, ps_stratum = fp3, mcmc = fp7, overall_rates = fp9, rates_by_stratum = fp10)
  }
  output <- list(n_s1 = n_s1, n_s0 = n_s0,
                 ne_s1 = ne_s1, ne_s0 = ne_s0,
                 r_s1 = r_s1, r_s0 = r_s0,
                 r_sat = R_SAT, r_eca = R_ECA,
                 v = v, r = r, A_eff = sum(A * r * gPPP),
                 PPP = PPP, gPPP = gPPP, gamma = gamma,
                 mean.theta = post.mean, var.theta = post.var, CrI = CrI, Width = Width,
                 theta_s = theta_s, theta_o = theta_o)
  if (!is.null(plots_list)) output$plots <- plots_list
  return(output)
}

