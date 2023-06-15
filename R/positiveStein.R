#' The positive rule Stein estimator
#'
#' This function calculates the positive-rule Stein estimatorm that is a improved version of Stein Estimator by only considering the positive part of shrinking factor
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
#' @param H A given q_restr x p_vars matrix.
#' @param h A given q_restr x 1 vector.
#'
#' @return A vector of regression coefficients
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
#'
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' positiveStein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' positiveStein(X, y, H, h)
#' @export
positiveStein <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  beta <- r_est + as.numeric(1 - d / test_stat) * as.integer(test_stat > d) * (u_est - r_est)
  residuals <- y-X%*%beta
  fit <- structure(list(coefficients=beta,residuals=residuals), class = c("positiveStein"))
  fit
}



#'Predict method for Linear Model Fits
#'
#'Predicted values based on linear model object.
#'
#' @param object An object of class "\code{positiveStein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{Stein}" or "\code{unrestricted}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#'
#' @seealso \code{\link{predict.Arima}}, \code{\link{predict.bats}},
#' \code{\link{predict.tbats}}, \code{\link{predict.ets}},
#' \code{\link{predict.nnetar}}
#' \code{\link{fitted.Arima}}, \code{\link{fitted.bats}},
#' \code{\link{fitted.tbats}}, \code{\link{fitted.ets}},
#' \code{\link{fitted.nnetar}}.
#'
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model<-positiveStein(X, y, H, h)
#' fitted(model)
#' @export
fitted.positiveStein <- function(object,newdata){
    return(newdata%*%object$beta)
}


residuals.positiveStein <- function(object){
    return(object$residuals)
}


coefficients.positiveStein <- function(object){
    return(object$coefficients)
}




