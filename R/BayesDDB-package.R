#' Bayesian Dynamic Borrowing with Double Adjustment
#'
#' @keywords internal
#' @name _PACKAGE
#' @rdname BayesDDB-package
#' @docType package
#' @importFrom grDevices dev.new
#' @importFrom ggplot2 geom_density geom_hline
#' @importFrom stats acf glm lag predict qbeta qnorm quantile rbeta rbinom sd update var
#' @description
#' BayesDDB implements a Bayesian dynamic borrowing approach for augmenting
#' clinical trial data with external data. The method uses propensity score
#' stratification for baseline adjustment and an elastic prior based on
#' posterior predictive probability for outcome adjustment.
#'
#' Main functions:
#' \itemize{
#'   \item \code{\link{bdb_single_arm_bin}} - Single-arm trial, binary outcome
#'   \item \code{\link{bdb_single_arm_cont}} - Single-arm trial, continuous outcome (requires JAGS)
#'   \item \code{\link{bdb_two_arms_bin}} - Two-arm trial, binary outcome
#'   \item \code{\link{bdb_two_arms_cont}} - Two-arm trial, continuous outcome (requires JAGS)
#'   \item \code{\link{elastic_function}} - Arctangent elastic function for outcome discounting
#' }
#'
#' See \code{vignette("Introduction to BayesDDB", package = "BayesDDB")} for methodology and examples.
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
"_PACKAGE"
