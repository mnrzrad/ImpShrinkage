#' The preliminary test estimator
#'
#' This function calculates the improved estimator after performing the preliminary test for the null hipothesis: H.beta = h
#'
#' @param X  Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y  Univariate quantitative response variable with dimension n_obs
#' @param H  A given q_restr x p_vars matrix.
#' @param h  A given q_restr x 1 vector.
#' @param alpha  A given significance level
#'
#' @return  A vector of regression coefficients
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' # H beta != h
#' p <- ncol(X)
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' @importFrom stats qf
#' @export

preliminaryTest <- function(X, y, H, h, alpha) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)

  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  threshold <- stats::qf(1 - alpha, q, n - p)
  u_est - (u_est - r_est) * as.integer(test_stat < threshold)
}
