#' The Unrestricted estimator
#'
#' This function calculates the unrestricted estimator
#'
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y Quantitative response variable.
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
#' unrestricted(X,y)
#'
#' @export
#'
unrestricted <- function(X, y) {
  solve(t(X) %*% X) %*% t(X) %*% y
}
