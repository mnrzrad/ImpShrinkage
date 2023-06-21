#' The restricted estimator
#'
#' This function calculates the restricted estimator using
#' \deqn{\hat{\beta}^{R} = \hat{\beta}^{U} - (X^{\top}X)^{-1}H^{\top}
#' (H(X^{\top}X)^{-1}H^{\top})^{-1}(H\hat{\beta}^{U}-h)}
#' where
#' \itemize{
#'   \item \eqn{\hat{\beta}^{U}} is the unrestricted estimator; See \code{\link{unrReg}}.
#'   \item \eqn{H\beta = h} represents a subspace of the parameter space induced
#' by the non-sample information. Here, \eqn{H} is a known \eqn{q \times p}
#' matrix, and \eqn{h} is a known \eqn{q}-vector.
#' }
#'
#'
#'#' The corresponding estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{R})^{\top}(y - X\hat{\beta}^{R}).}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Vector with response observations of size \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#'
#'
#' @returns
#' An object of class \code{restricted} is a list containing at least the following components:
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
#' Kaciranlar, S., Akdeniz, S. S. F., Styan, G. P. & Werner, H. J. (1999). A new biased
#' estimators in linear regression and detailed
#' analysis of the widely-analysed dataset on
#' portland cement. \emph{Sankhya, Series B}, 61(3), 443-459.
#'
#' Kibria, B. M. Golam (2005). Applications of Some Improved Estimators in Linear Regression,
#' \emph{Journal of Modern Applied Statistical Methods}, 5(2), 367- 380.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' resReg(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' resReg(X, y, H, h)
#'
#'
#' data(cement)
#' X <- as.matrix(cbind(1, cement[, 1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0, 1, -1, 1, 0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' resReg(X, y, H, h)
#' # Based on Kibria (2005)
#' H <- matrix(c(0, 1, -1, 1, 0, 0, 0, 1, -1, -1, 0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' resReg(X, y, H, h)
#' @export

resReg <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  u_est <- unrReg(X, y) # unrestricted estimator
  C <- t(X) %*% X
  beta <- u_est$coef - solve(C) %*% t(H) %*% solve(H %*% solve(C) %*% t(H)) %*% (H %*% u_est$coef - h)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, s2 = s2, residuals = residuals, fitted.value = fittedValues), class = c("restricted"))
  fit
}

#' Extract Model Fitted Values
#'
#' Fitted values based on object \code{restrcited}.
#'
#' @param object An object of class \code{restricted}.
#' @param ... Other arguments.
#'
#' @return Fitted values extracted from the object \code{restricted}.
#'
#' @seealso
#' \code{\link{fitted.unrestricted}},
#' \code{\link{fitted.preliminaryTest}},
#' \code{\link{fitted.improvedpreliminaryTest}},
#' \code{\link{fitted.stein}},
#' \code{\link{fitted.positivestein}}
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- resReg(X, y, H, h)
#' fitted(model)
#' @export
fitted.restricted <- function(object, ...) {
  return(object$fitted.value)
}

#' Extract Model Predictions Values
#'
#' Predicted values based on object \code{restrcited}.
#'
#' @param object An object of class \code{restricted}.
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other arguments.
#'
#' @return A vector of predictions.
#'
#' @seealso
#' \code{\link{predict.unrestricted}},
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- resReg(X, y, H, h)
#' predict(model, X)
#' @export
predict.restricted <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' Extract Model Residuals
#'
#' Residuals values based on model object \code{restricted}.
#'
#' @param object An object of class \code{restricted}.
#' @param ... Other arguments.
#'
#' @return A vector of residuals.
#'
#' \code{\link{residuals.unrestricted}},
#' \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.improvedpreliminaryTest}},
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- resReg(X, y, H, h)
#' residuals(model)
#' @export

residuals.restricted <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' Coefficients extracted from the model object \code{restrcited}.
#'
#' @param object An object of class \code{restricted}.
#' @param ... Other arguments.
#'
#' @return A vector of coefficients.
#'
#' @seealso
#' \code{\link{coefficients.unrestricted}},
#' \code{\link{coefficients.preliminaryTest}},
#' \code{\link{coefficients.improvedpreliminaryTest}},
#' \code{\link{coefficients.stein}},
#' \code{\link{coefficients.positivestein}},
#' \code{\link{coef.unrestricted}},
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- resReg(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.restricted <- function(object, ...) {
  return(object$coef)
}

#' @rdname coefficients.restricted
#' @importFrom stats coef
#' @examples
#' coef(model)
#' @export

coef.restricted <- function(object, ...) {
  return(object$coef)
}
