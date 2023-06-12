#' Test-Statistics
#'
#' This function calculates the test statistics
#'
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y Quantitative response variable.
#' @param H A given q x p matrix.
#' @param h A given q x 1 vector.
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
#' h <- rep(0, nrow(H))
#' test_statistics(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' test_statistics(X, y, H, h)
#' @export
test_statistics <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- n - p
  q <- nrow(H)
  u_est <- unrestricted(X, y) # unrestricted estimator
  s2 <- (t(y - X %*% u_est) %*% (y - X %*% u_est)) / m
  C <- t(X) %*% X
  diff <- (H %*% u_est - h)
  as.numeric((t(diff) %*% solve(H %*% solve(C) %*% t(H)) %*% diff) / (q * s2))
}
