#' Arctangent Elastic Function
#'
#' Applies the arctangent elastic function to transform the posterior predictive
#' probability (PPP). This function is used for outcome adjustment in the
#' double-adjustment approach.
#'
#' @param ppp Numeric vector of posterior predictive probabilities (values
#'        between 0 and 1).
#' @param a Numeric. Tuning parameter for the arctangent elastic function
#'        (must be >= 1). Higher values result in less discounting when PPP
#'        is near 0.5 and more aggressive discounting when PPP deviates from 0.5.
#'
#' @return Numeric vector of transformed probabilities, same length as \code{ppp}.
#'
#' @details
#' The elastic function is defined as:
#' \deqn{g(PPP) = \frac{\arctan(a \cdot \sin(\pi \cdot PPP))}{\arctan(a)}}
#'
#' This transformation has the following properties:
#' \itemize{
#'   \item When PPP = 0.5, g(PPP) = 1 (no discounting - outcomes are compatible)
#'   \item When PPP approaches 0 or 1, g(PPP) approaches 0 (maximum discounting)
#'   \item The function is symmetric around PPP = 0.5
#'   \item Parameter \code{a} controls the steepness of discounting
#' }
#'
#' @examples
#' # Example with different PPP values
#' ppp_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' elastic_function(ppp_values, a = 1)
#' elastic_function(ppp_values, a = 5)
#'
#' # Plot the elastic function for different tuning parameters
#' x <- seq(0, 1, by = 0.01)
#' plot(x, elastic_function(x, a = 1), type = "l",
#'      xlab = "PPP", ylab = "g(PPP)", main = "Elastic Function")
#' lines(x, elastic_function(x, a = 5), col = "red")
#' lines(x, elastic_function(x, a = 10), col = "blue")
#' legend("bottomright", legend = c("a=1", "a=5", "a=10"),
#'        col = c("black", "red", "blue"), lty = 1)
#'
#' @export
elastic_function <- function(ppp, a = 1) {
  if (!is.numeric(ppp)) {
    stop("'ppp' must be numeric")
  }
  if (any(ppp < 0 | ppp > 1, na.rm = TRUE)) {
    stop("'ppp' values must be between 0 and 1")
  }
  if (!is.numeric(a) || length(a) != 1 || a < 1) {
    stop("'a' must be a single numeric value >= 1")
  }

  atan(a * sin(ppp * pi)) / atan(a)
}


#' Calculate Mean and Confidence Interval
#'
#' Computes the sample mean and confidence interval for a numeric vector.
#'
#' @param x Numeric vector of observations.
#' @param conf_level Numeric. Confidence level for the interval (default 0.95).
#'
#' @return Named numeric vector with three elements:
#'         \code{mean}, \code{lower}, and \code{upper}.
#'
#' @details
#' The confidence interval is computed using the normal approximation:
#' \deqn{\bar{x} \pm z_{\alpha/2} \cdot \frac{s}{\sqrt{n}}}
#'
#' @examples
#' x <- rnorm(100, mean = 5, sd = 2)
#' mean_ci(x)
#' mean_ci(x, conf_level = 0.99)
#'
#' @keywords internal
#' @noRd
mean_ci <- function(x, conf_level = 0.95) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric")
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1")
  }

  x <- x[!is.na(x)]
  n <- length(x)

  if (n == 0) {
    return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }

  mu <- mean(x)
  sigma <- sd(x)
  za <- qnorm(1 - (1 - conf_level) / 2)
  se <- sigma / sqrt(n)

  c(mean = mu, lower = mu - za * se, upper = mu + za * se)
}


#' Mean and CI by significance level (internal)
#' @param x numeric vector
#' @param alpha significance level (e.g. 0.05 for 95\% CI)
#' @return 3-vector (mean, lower, upper)
#' @keywords internal
#' @noRd
mean_CI <- function(x, alpha) {
  res <- mean_ci(x, conf_level = 1 - alpha)
  c(res[["mean"]], res[["lower"]], res[["upper"]])
}
