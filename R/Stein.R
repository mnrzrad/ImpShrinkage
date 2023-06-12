#' The Stein estimator
#'
#' This function calculates the Stein estimator
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
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' Stein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' Stein(X, y, H, h)
#' @export

Stein <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  u_est - d * (u_est - r_est) / test_stat
}

