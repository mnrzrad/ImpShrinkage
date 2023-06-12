#' The preliminary test estimator
#'
#' This function calculates the preliminary test estimator
#'
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y Qunatitavie response variable.
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
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, q)
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' @export

preliminaryTest <- function(X, y, H, h, alpha) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)

  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h, q)
  threshold <- qf(1 - alpha, q, n - p)
  u_est - (u_est - r_est) * as.integer(test_stat < threshold)
}


