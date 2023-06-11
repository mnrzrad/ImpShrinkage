#' The restricted estimator
#'
#' positiveStein
#'
#' @param x scaler.
#' @param q scaler.
#' @param n number.
#' @param alpha float.
#'
#' @references ref
#' @examples
#' n <- 100
#' p <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simulate_data(n, p, beta)
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
