#' Test-Statistics
#'
#' This function calculates the test statistics used in the "\code{preliminaryTest}".
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
#' @param H A given q_restr x p_vars matrix.
#' @param h A given q_restr x 1 vector.
#'
#' @return A numerical value of the test statistic
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
#'
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' test_statistics(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
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
