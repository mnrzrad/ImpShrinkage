#' The Unrestricted estimator
#'
#' This function calculates the unrestricted estimator that may also be used in the preliminaryTest.R
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
#'
#' @return A vector of regression coefficients
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
#'
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' unrestricted(X, y)
#'
#' @export
#'
unrestricted <- function(X, y) {
  solve(t(X) %*% X) %*% t(X) %*% y
}
