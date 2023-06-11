#' The preliminary test estimator
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


