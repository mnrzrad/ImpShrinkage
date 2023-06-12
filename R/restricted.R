#' The restricted estimator
#'
#' This function calculates the restricted estimator
#'
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y Quantitative response variable.
#' @param H A given q x p matrix.
#' @param h A given q x 1 vector.
#'
#' @return A vector of coefficients
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
#'
#' @examples
#' n <- 100
#' p <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n, p, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, q)
#' restricted(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, q)
#' restricted(X, y, H, h)
#'
#' @export

restricted <- function(X, y, H, h) {
  u_est <- unrestricted(X, y) # unrestricted estimator
  C <- t(X) %*% X
  u_est - solve(C) %*% t(H) %*% solve(H %*% solve(C) %*% t(H)) %*% (H %*% u_est - h)
}
