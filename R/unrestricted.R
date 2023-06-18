#' The Unrestricted estimator
#'
#' This function calculates the unrestricted estimator using the formula:
#' \deqn{\hat{\beta}^{U} = (X^{\top} X)^{-1} X^{\top} y}
#' where \eqn{^{\top}} denotes the transpose of a matrix, calculated using
#' \code{\link[base]{t}}, and \eqn{^{-1}} represents the inverse of a matrix,
#' calculated using \code{\link[base]{solve}}. It is important to note that the
#' input matrices \eqn{X} and \eqn{y} should be standardized, for example, by
#' using \code{\link[base]{scale}}. Alternatively, the user can employ
#' \code{\link[stats]{lm}} to obtain this estimator, but it is crucial to
#' remember to set \code{intercept = FALSE}.
#'
#' The corresponding unrestricted estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{U})^{\top}(y - X\hat{\beta}^{U})}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#'
#' @param y Univariate quantitative response variable with dimension \code{n}.
#'
#' @return A vector of regression coefficients
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and
#'  Stein‐Type Estimation With Applications}, Wiley.
#'
#'  ref
#'
#' @examples
#' data(cement)
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' unrestricted(X, y) #'
#' @export
#'
unrestricted <- function(X, y) {
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  residuals <- (y - X %*% beta)[, 1]
  n <- dim(X)[1]
  p <- dim(X)[2]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(
    coef = beta, residuals = residuals, s2 = s2,
    fitted.value = fittedValues
  ), class = c("unrestricted"))
  fit
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
#' fitted(model)
#' @export
fitted.unrestricted <- function(object, ...) {
  return(object$fitted.value)
}

#' @rdname predict.stein
#' @importFrom stats predict
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
