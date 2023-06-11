#' Simulation data
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
#'
#' @export

simulate_data <- function(n, p, beta){
    X <- matrix(rnorm(n*p), nr = n, nc = p)
    y <- X %*% beta + rnorm(n)
    return(list(X= X, y = y))
}
