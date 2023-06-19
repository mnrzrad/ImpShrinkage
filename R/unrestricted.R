#' The Unrestricted estimator
#'
#' This function calculates the unrestricted estimator as
#' \deqn{\hat{\beta}^{U} = (X^{\top} X)^{-1} X^{\top} y}
#' where \eqn{^{\top}} denotes the transpose of a matrix. It is important to note that the
#' input matrices \eqn{X} and \eqn{y} should be standardized, for example, by
#' using \code{\link[base]{scale}}. Alternatively, the user can employ
#' \code{\link[stats]{lm}} to obtain this estimator, but it is crucial to
#' remember to set \code{intercept = FALSE}.
#'
#' The corresponding unrestricted estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{U})^{\top}(y - X\hat{\beta}^{U})}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p}, where
#' each row is an observation vector.
#' @param y Vector with response observations of size \code{n}.
#'
#' @returns
#' An object of class \code{unrestricted} is a list containing at least the following components:
#'   \describe{
#'     \item{\code{coef}}{A named vector of coefficients.}
#'     \item{\code{residuals}}{The residuals, that is, response minus fitted values.}
#'     \item{\code{s2}}{The estimated variance.}
#'     \item{\code{fitted.values}}{The fitted values.}
#'   }
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and
#'  Stein‐Type Estimation With Applications}, Wiley.
#'
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

#' Extract Model Fitted Values
#'
#' \code{fitted} is a generic function which extracts fitted values from objects
#'  returned by modeling functions. \code{fitted.values} is an alias for it.
#'
#' @param object An object of class "\code{unrestricted}".
#' @param ... Other.
#' @seealso#'
#' \code{\link{fitted.restricted}},
#' \code{\link{fitted.preliminaryTest}},
#' \code{\link{fitted.improvedpreliminaryTest}},
#' \code{\link{fitted.stein}},
#' \code{\link{fitted.positivestein}}.
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
#'
fitted.unrestricted <- function(object, ...) {
  return(object$fitted.value)
}

#' Model Predictions
#'
#' \code{predict} is a generic function for predictions from the results of various
#' model fitting functions.
#'
#' @param object An object of class "\code{unrestricted}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other.
#' @seealso
#' \code{\link{predict.restricted}},
#' \code{\link{predict.preliminaryTest}},
#' \code{\link{predict.improvedpreliminaryTest}},
#' \code{\link{predict.stein}},
#' \code{\link{predict.positivestein}}.
#'
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
#'
predict.unrestricted <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' residuals method for Model Fits
#'
#' residuals values based on model object.
#'
#' @param object An object of class "\code{unrestricted}".
#' @param ... Other.
#' @seealso
#' \code{\link{residuals.restricted}},
#' \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.improvedpreliminaryTest}}
#' \code{\link{residuals.stein}},
#' \code{\link{residuals.positivestein}}.
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

#' Extract Model Coefficients
#'
#' \code{coef} is a generic function which extracts model
#' coefficients from objects returned by modeling \code{functions.coefficients}
#' is an alias for it.
#'
#' @param object An object of class "\code{unrestricted}".
#' @param ... Other.
#' @seealso
#' \code{\link{coefficients.restricted}},
#' \code{\link{coefficients.preliminaryTest}},
#' \code{\link{coefficients.improvedpreliminaryTest}},
#' \code{\link{coefficients.stein}},
#' \code{\link{coefficients.positivestein}},
#' \code{\link{coef.restricted}},
#' \code{\link{coef.preliminaryTest}},
#' \code{\link{coef.improvedpreliminaryTest}}
#' \code{\link{coef.stein}},
#' \code{\link{coef.positivestein}}.
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

#' @rdname coefficients.unrestricted
#' @importFrom stats coef
#' @export
#' @examples
#' coefficients(model)
coef.unrestricted <- function(object, ...) {
  return(object$coef)
}
