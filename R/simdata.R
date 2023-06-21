#' Simulation data
#'
#' This function generates a toy example. The error term, \eqn{\varepsilon},
#' and the design matrix, \eqn{X}, are simulated from standard normal
#' distributions, \eqn{\mathcal{N}(0,1)}, using the \code{\link[stats]{rnorm}}
#' function. Given the true parameter vector, \eqn{\beta}, the response vector,
#' \eqn{y}, is calculated as
#' \deqn{y = X \beta + \varepsilon.}
#'
#'
#' @param n Number of observations.
#'
#' @param p Number of variables.
#'
#' @param beta Regression parameter.
#'
#' @param seed (Optional) The random seed for reproducibility. Default is \code{NULL}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{X}{a matrix of dimensions \code{n x p}.}
#'   \item{y}{a numeric vector of length \code{n}.}
#' }
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and
#'  Stein‐Type Estimation With Applications}, Wiley.
#'
#'
#' @examples
#' simulated_data <- simdata(n = 100, p = 5, beta = c(2, 1, 3, 0, 5))
#' X <- simulated_data$X
#' y <- simulated_data$y
#' X
#' y
#'
#' @importFrom stats rnorm
#' @export

simdata <- function(n, p, beta, seed = NULL) {
  set.seed(seed)

  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)

  y <- X %*% beta + stats::rnorm(n)

  return(list(X = X, y = y))
}
