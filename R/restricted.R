#' The restricted estimator
#'
#' This function calculates the restricted estimator that may also be used in the preliminaryTest.R
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
#' @param H A given q_restr x p_vars matrix.
#' @param h A given q_restr x 1 vector.
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
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' restricted(X, y, H, h)
#'
#' # H beta != h
#' p <- ncol(X)
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' restricted(X, y, H, h)
#'
#' @export

restricted <- function(X, y, H, h) {
  u_est <- unrestricted(X, y) # unrestricted estimator
  C <- t(X) %*% X
  u_est - solve(C) %*% t(H) %*% solve(H %*% solve(C) %*% t(H)) %*% (H %*% u_est - h)
}
