#' Simulation data
#'
#' This function generates simulated data from a standard normal distribution as an example.
#'
#' @param n Number of observations.
#'   An integer value specifying the total number of observations in the dataset.
#'   The function expects a positive integer greater than zero.
#'
#' @param p Number of variables.
#'   An integer value indicating the number of predictor variables in the regression model.
#'   The function expects a positive integer greater than zero and more than \code{n}.
#'
#' @param beta Regression parameter.
#'   A numeric vector of length p representing the regression coefficients for the predictor variables.
#'   Each element in the vector corresponds to the coefficient for a specific predictor variable.
#'   The function uses these coefficients to estimate the relationship between the predictor variables and the response variable.
#'   The length of the vector must be equal to the value specified in the \code{p} parameter.
#'   If not provided, default values may be used.
#'
#' @param seed (Optional) The random seed for reproducibility.
#'   An integer value representing the seed for the random number generator. Default is NULL.
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
#' simulated_data <- simulate(n = 100, p = 5)
#' X <- simulated_data$X
#' y <- simulated_data$y
#'
#' simulated_data <- simulate(n = 100, p = 5, beta = c(2, 1, 3, 0, 5))
#' X <- simulated_data$X
#' y <- simulated_data$y
#'
#' @export

simdata <- function(n, p, beta = rep(1, p), seed = NULL) {
  set.seed(seed)

  X <- matrix(rnorm(n * p), nr = n, nc = p)

  y <- X %*% beta + rnorm(n)

  return(list(X = X, y = y))
}
