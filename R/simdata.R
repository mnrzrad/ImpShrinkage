#' Simulation data
#'
#' This function generates simulated data from a standard normal distribution as an example.
#'
#' @param n Number of observations.
#'
#' @param p Number of variables.
#'
#' @param beta Regression parameter.
#'
#'
#' @param seed (Optional) The random seed for reproducibility. Default is NULL.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{X}{a matrix of dimensions \code{n x p}}
#'   \item{y}{a numeric vector of length \code{n}}
#' }
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
#'
#' @examples
#' simulated_data <- simdata(n = 100, p = 5)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' X
#' y
#' simulated_data <- simdata(n = 100, p = 5, beta = c(2, 1, 3, 0, 5))
#' X <- simulated_data$X
#' y <- simulated_data$y
#' X
#' y
#' @export

simdata <- function(n, p, beta, seed = NULL) {
  set.seed(seed)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)

  y <- X %*% beta + rnorm(n)

  return(list(X = X, y = y))
}
