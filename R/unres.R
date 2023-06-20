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
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{U})^{\top}(y - X\hat{\beta}^{U}).}
#'
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p}, where
#' each row is an observation vector;
#' @param y Vector with response observations of size \code{n}.
#'
#' @returns
#' An object of class \code{unrestricted} is a list containing at least the following components:
#'   \describe{
#'     \item{\code{coef}}{A named vector of coefficients.}
#'     \item{\code{residuals}}{The residuals, that is, the response values minus fitted values.}
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
#' unres(X, y)
#'
#' data(cement)
#' X <- as.matrix(cbind(1, cement[, 1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0, 1, -1, 1, 0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' unres(X, y)
#'
#' H <- matrix(c(0, 1, -1, 1, 0, 0, 0, 1, -1, -1, 0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' unres(X, y)
#' @export
#'
unres <- function(X, y) {
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  residuals <- (y - X %*% beta)[, 1]
  n <- dim(X)[1]
  p <- dim(X)[2]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, s2 = s2, residuals = residuals, fitted.value = fittedValues), class = c("unrestricted"))
  fit
}

#' Extract Model Fitted Values
#'
#' Fitted values based on object \code{unrestrcited}.
#'
#' @param object An object of class \code{unrestricted}.
#' @param ... Other arguments.
#'
#' @return A vector of fitted values.
#'
#' @seealso
#' \code{\link{fitted.res}},
#' \code{\link{fitted.pt}},
#' \code{\link{fitted.ipt}},
#' \code{\link{fitted.st}},
#' \code{\link{fitted.pst}}.
#'
#' @importFrom stats fitted
#'
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

#' Extract Model Predictions Values
#'
#' Predicted values based on object \code{unrestrcited}.
#'
#' @param object An object of class \code{unrestricted}.
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other arguments.
#'
#' @return A vector of predictions.
#'
#' @seealso
#' \code{\link{predict.res}},
#' \code{\link{predict.pt}},
#' \code{\link{predict.ipt}},
#' \code{\link{predict.st}},
#' \code{\link{predict.pst}}.
#'
#' @importFrom stats predict
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unres(X, y)
#' predict(model, X)
#' @export
#'
predict.unres <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' Extract Model Residuals
#'
#' Residuals values based on model object \code{unrestrcited}.
#'
#' @param object An object of class \code{unrestricted}.
#' @param ... Other arguments.
#'
#' @retrun A vector of residuals.
#'
#' @seealso
#' \code{\link{residuals.res}},
#' \code{\link{residuals.pt}},
#' \code{\link{residuals.ipt}}
#' \code{\link{residuals.st}},
#' \code{\link{residuals.pst}}.
#' @importFrom stats residuals
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unres(X, y)
#' residuals(model)
#' @export

residuals.unres <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' Coefficients extracted from the model object \code{unrestrcited}.
#'
#' @param object An object of class \code{unrestricted}.
#' @param ... Other arguments.
#'
#' @return A vector of coefficients.
#'
#' @seealso
#' \code{\link{coefficients.res}},
#' \code{\link{coefficients.pt}},
#' \code{\link{coefficients.ipt}},
#' \code{\link{coefficients.st}},
#' \code{\link{coefficients.pst}},
#' \code{\link{coef.res}},
#' \code{\link{coef.pt}},
#' \code{\link{coef.ipt}}
#' \code{\link{coef.st}},
#' \code{\link{coef.pst}}.
#' @importFrom stats coefficients
#' @examples
#' n_obs <- 100
#' p_vars <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simdata(n = n_obs, p = p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' model <- unres(X, y)
#' coefficients(model)
#' @export

coefficients.unres <- function(object, ...) {
  return(object$coef)
}

#' @rdname coefficients.unres
#' @importFrom stats coef
#' @export
#' @examples
#' coef(model)
coef.unres <- function(object, ...) {
  return(object$coef)
}
