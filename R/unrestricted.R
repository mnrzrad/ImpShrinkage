#' The Unrestricted estimator
#'
#' This function calculates the unrestricted estimator that may also be used in the "\code{\link{preliminaryTest}}".
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
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
#' unrestricted(X, y)
#'
#' @export
#'
unrestricted <- function(X, y) {
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  residuals <- (y - X %*% beta)[, 1]
  d <- dim(X)
  s2 <- sum(residuals^2) / (d[1] - d[2])
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2), class = c("unrestricted"))
  fit
}


#' @rdname fitted.stein
#' @importFrom stats predict
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unrestricted(X, y)
#' fitted(model, X)
#' @export
fitted.unrestricted <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' @rdname fitted.stein
#' @importFrom stats fitted
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unrestricted(X, y)
#' predict(model, X)
#' @export
predict.unrestricted <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' @rdname residuals.stein
#' @importFrom stats residuals
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unrestricted(X, y)
#' residuals(model)
#' @export

residuals.unrestricted <- function(object, ...) {
  return(object$residuals)
}

#' @rdname coefficients.stein
#' @importFrom stats coefficients
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unrestricted(X, y)
#' coefficients(model)
#' @export

coefficients.unrestricted <- function(object, ...) {
  return(object$coef)
}

#' @rdname coefficients.stein
#' @importFrom stats coef
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unrestricted(X, y)
#' coef(model)
#' @export

coef.unrestricted <- function(object, ...) {
  return(object$coef)
}
