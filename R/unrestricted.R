#' The Unrestricted estimator
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
#' unrestricted(X,y)
#'
#' @export
#'
unrestricted <- function(X, y) {
  solve(t(X) %*% X) %*% t(X) %*% y
}
