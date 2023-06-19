#' The positive-rule Stein estimator
#'
#' This function calculates the positive-rule Stein estimator. This estimator is an improved
#' version of the Stein estimator, where only the positive part of the
#' shrinking factor is considered. It may be calculated by
#' \deqn{\hat{\beta}^{S+}= \hat{\beta}^{S} + (1 + d \mathcal{L}^{-1}) I(\mathcal{L} > d) (\hat{\beta}^{U} - \hat{\beta}^{R})}
#' where \eqn{I(A)} denotes an indicator function and
#' \itemize{
#'   \item \eqn{\hat{\beta}^{S}} is the \code{\link{stein}} estimator;
#'   \item \eqn{\hat{\beta}^{U}} is the \code{\link{unrestricted}} estimator;
#'   \item \eqn{\hat{\beta}^{R}} is the \code{\link{restricted}} estimator;
#'   \item \eqn{\mathcal{L}} is the \code{\link{test_statistic}};
#'   \item and \eqn{d} is the shrinkage factor.
#' }
#'
#' The corresponding estimator of \eqn{\sigma^2} is given by
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{S+})^{\top}(y - X\hat{\beta}^{S+})}
#'
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Vector with response observations of size \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#' @param d (optional) If not provided (or set to \code{NULL}), it will be
#' calculated using \eqn{\frac{{(q - 2) \cdot (n - p)}}{{q \cdot (n - p + 2)}}}
#' @param is_error_normal logical value indicating whether the errors follow a
#' normal distribution. If \code{is_error_normal} is \code{TRUE}, the distribution
#' of the test statistics for the null hypothesis is F distribution
#' (\code{\link[stats]{FDist}}). On the other hand, if the errors have a
#' non-normal distribution, the asymptotic distribution of the test statistics
#' is \eqn{\chi^2} distribution (\code{\link[stats]{Chisquare}}). By default,
#' \code{is_error_normal} is set to \code{FALSE}.
#'
#'
#' @returns
#' An object of class \code{positivestein} is a list containing at least the following components:
#'   \describe{
#'     \item{\code{coef}}{A named vector of coefficients.}
#'     \item{\code{residuals}}{The residuals, that is, the response values minus fitted values.}
#'     \item{\code{s2}}{The estimated variance.}
#'     \item{\code{fitted.values}}{The fitted values.}
#'   }
#'
#'
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test
#'   and Stein‐Type Estimation With Applications}, Wiley.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' positivestein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' positivestein(X, y, H, h)
#'
#' data(cement)
#' X <- as.matrix(cbind(1, cement[, 1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0, 1, -1, 1, 0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' positivestein(X, y, H, h)
#'
#' H <- matrix(c(0, 1, -1, 1, 0, 0, 0, 1, -1, -1, 0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' positivestein(X, y, H, h)
#'
#' @export
positivestein <- function(X, y, H, h, d = NULL, is_error_normal = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistic(X, y, H, h, is_error_normal = is_error_normal)
  beta <- r_est$coef + as.numeric(1 - d / test_stat) * as.integer(test_stat > d) * (u_est$coef - r_est$coef)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2, fitted.value = fittedValues), class = c("positivestein"))
  fit
}


#' Extract Model Fitted Values
#'
#' \code{fitted} is a generic function which extracts fitted values from objects
#'  returned by modeling functions. \code{fitted.values} is an alias for it.
#'
#' @param object An object of class "\code{positivestein}".
#' @param ... Other.
#' @seealso#' \code{\link{fitted.unrestricted}}, \code{\link{fitted.restricted}},
#' \code{\link{fitted.preliminaryTest}},\code{\link{fitted.improvedpreliminaryTest}},
#' \code{\link{fitted.stein}}, \code{\link{fitted.positivestein}}
#' @importFrom stats fitted
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
#' model <- positivestein(X, y, H, h)
#' fitted(model)
#' @export
fitted.positivestein <- function(object, ...) {
  return(object$fitted.value)
}

#' Model Predictions
#'
#' \code{predict} is a generic function for predictions from the results of various
#' model fitting functions.
#'
#' @param object An object of class \code{positivestein}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other.
#' @seealso \code{\link{predict.unrestricted}}, \code{\link{predict.restricted}},
#'  \code{\link{predict.preliminaryTest}}, \code{\link{predict.improvedpreliminaryTest}},
#'  \code{\link{predict.stein}}, \code{\link{predict.positivestein}}.
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
#' model <- positivestein(X, y, H, h)
#' predict(model, X)
#' @export
predict.positivestein <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' residuals method for Model Fits
#'
#' residuals values based on model object.
#'
#' @param object An object of class "\code{positivestein}"
#' @param ... Other.
#' @seealso \code{\link{residuals.positivestein}}, \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.restricted}}, \code{\link{residuals.stein}},
#' \code{\link{residuals.unrestricted}}, \code{\link{residuals.improvedpreliminaryTest}}.
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
#' model <- positivestein(X, y, H, h)
#' residuals(model)
#' @export

residuals.positivestein <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' \code{coef} is a generic function which extracts model
#' coefficients from objects returned by modeling \code{functions.coefficients}
#' is an alias for it.
#'
#' @param object An object of class "\code{positivestein}".
#' @param ... Other.
#' @seealso
#' \code{\link{coefficients.unrestricted}},
#' \code{\link{coefficients.restricted}},
#' \code{\link{coefficients.preliminaryTest}},
#' \code{\link{coefficients.improvedpreliminaryTest}},
#' \code{\link{coefficients.stein}},
#' \code{\link{coef.unrestricted}},
#' \code{\link{coef.restricted}},
#' \code{\link{coef.preliminaryTest}},
#' \code{\link{coef.improvedpreliminaryTest}}
#' \code{\link{coef.stein}}.
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
#' model <- positivestein(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.positivestein <- function(object, ...) {
  return(object$coef)
}

#' @rdname coefficients.positivestein
#' @importFrom stats coef
#' @examples
#' coef(model)
#' @export

coef.positivestein <- function(object, ...) {
  return(object$coef)
}
