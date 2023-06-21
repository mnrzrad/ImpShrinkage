#' The improved preliminary test estimator
#'
#' This function calculates the improved preliminary test estimator. When the error has a normal distribution,
#' this estimator can be calculated by
#' \deqn{\hat{\beta}^{iPT}= \hat{\beta}^{PT} - d (\hat{\beta}^{U} - \hat{\beta}^{R}) \mathcal{L}^{-1} I(\mathcal{L} > F_{q,n-p}(\alpha)) }
#' and, when the error has a non-normal distribution, by
#' \deqn{\hat{\beta}^{iPT}= \hat{\beta}^{PT} - d (\hat{\beta}^{U} - \hat{\beta}^{R}) \mathcal{L}^{-1} I(\mathcal{L} > \chi^2_{q}(\alpha)),}
#' where \eqn{I(A)} denotes an indicator function and
#' \itemize{
#'   \item \eqn{\hat{\beta}^{PT}} is the preliminary test estimator; See \link[preliminaryTest]{pt}
#'   \item \eqn{\hat{\beta}^{U}} is the unrestricted estimator; See \link[unrestricted]{unres}.
#'   \item \eqn{\hat{\beta}^{R}} is the restricted estimator; See \link[restricted]{res}.
#'   \item \eqn{\mathcal{L}} is the test statistic; See \link{ImpShrinkage}[teststat].
#'   \item \eqn{F_{q,n-p}(\alpha)} is the upper \eqn{\alpha} level critical value of \eqn{F}-distribution with \eqn{(q,n-p)} degrees of freedom, calculated using \code{\link[stats]{qf}};
#'   \item \eqn{\chi^2_{q}(\alpha)} is the upper \eqn{\alpha} level critical value of \eqn{\chi^2}-distribution with \eqn{q} degree of freedom, calculated using \code{\link[stats]{qchisq}};
#'   \item \eqn{d} is the shrinkage factor;
#'   \item \eqn{\alpha} is the significance level.
#' }
#'
#' The corresponding estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{iPT})^{\top}(y - X\hat{\beta}^{iPT}).}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Vector with response observations of size \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#' @param alpha  A given significance level.
#' @param d (optional) If not provided (or set to \code{NULL}), it will be
#' calculated using \eqn{\frac{{(q - 2) \cdot (n - p)}}{{q \cdot (n - p + 2)}}.}
#' @param is_error_normal logical value indicating whether the errors follow a
#' normal distribution. If \code{is_error_normal} is \code{TRUE}, the distribution
#' of the test statistics for the null hypothesis is F distribution
#' (\code{\link[stats]{FDist}}). On the other hand, if the errors have a
#' non-normal distribution, the asymptotic distribution of the test statistics
#' is \eqn{\chi^2} distribution (\code{\link[stats]{Chisquare}}). By default,
#' \code{is_error_normal} is set to \code{FALSE}.
#'
#' @returns
#' An object of class \code{improvedpreliminaryTest} is a list containing at least the following components:
#'   \describe{
#'     \item{\code{coef}}{A named vector of coefficients.}
#'     \item{\code{residuals}}{The residuals, that is, the response values minus fitted values.}
#'     \item{\code{s2}}{The estimated variance.}
#'     \item{\code{fitted.values}}{The fitted values.}
#'   }
#'
#' @references
#' Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test and Stein‐Type Estimation With Applications}, Wiley.
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
#' ipt(X, y, H, h, alpha = 0.05)
#'
#' # H beta != h
#' p <- ncol(X)
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' ipt(X, y, H, h, alpha = 0.05)
#'
#' data(cement)
#' X <- as.matrix(cbind(1, cement[, 1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0, 1, -1, 1, 0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' ipt(X, y, H, h, alpha = 0.05)
#' # Based on Kibria (2005)
#' H <- matrix(c(0, 1, -1, 1, 0, 0, 0, 1, -1, -1, 0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' ipt(X, y, H, h, alpha = 0.05)
#' @importFrom stats qf
#' @export

ipt <- function(X, y, H, h, alpha, d = NULL, is_error_normal = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  q <- nrow(H)
  m <- n - p
  if (is.null(d)) {
    d <- ((q - 2) * m) / (q * (m + 2))
  }
  u_est <- unres(X, y)
  r_est <- res(X, y, H, h)
  pt_est <- pt(X, y, H, h, alpha)
  test_stat <- teststat(X, y, H, h, is_error_normal = is_error_normal)
  if (!is_error_normal) {
    threshold <- stats::qf(1 - alpha, q, n - p)
  } else {
    threshold <- stats::qchisq(1 - alpha, q)
  }
  threshold <- stats::qf(1 - alpha, q, n - p)
  beta <- pt_est$coef - d * ((u_est$coef - r_est$coef) / test_stat) * as.integer(test_stat < threshold)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, s2 = s2, residuals = residuals, fitted.value = fittedValues), class = c("improvedpreliminaryTest"))
  fit
}

#' Extract Model Fitted Values
#'
#' Fitted values based on object \code{improvedpreliminaryTest}.
#'
#' @param object An object of class \code{improvedpreliminaryTest}.
#' @param ... Other arguments.
#'
#' @return A vector of fitted values.
#'
#' @seealso
#' \code{\link{fitted.unres}},
#' \code{\link{fitted.res}},
#' \code{\link{fitted.pt}},
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- ipt(X, y, H, h, alpha = 0.05)
#' fitted(model)
#'
#' @export
fitted.ipt <- function(object, ...) {
  return(object$fitted.value)
}

#' Extract Model Predictions Values
#'
#' Predicted values based on object \code{improvedpreliminaryTest}.
#'
#' @param object An object of class "\code{improvedpreliminaryTest}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other arguments.
#'
#' @return A vector of predictions.
#'
#' @seealso
#' \code{\link{predict.unres}},
#' \code{\link{predict.res}},
#' \code{\link{predict.pt}},
#' \code{\link{predict.st}},
#' \code{\link{predict.pst}}.
#'
#' @importFrom stats predict
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
#' model <- ipt(X, y, H, h, alpha = 0.05)
#' predict(model, X)
#'
#' @export
#'
predict.ipt <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' Extract Model Residuals
#'
#' Residuals values based on model object \code{improvedpreliminaryTest}.
#'
#' @param object An object of class \code{improvedpreliminaryTest}.
#' @param ... Other arguments.
#'
#' @return A vector of residuals.
#'
#' @seealso
#' \code{\link{residuals.unres}},
#' \code{\link{residuals.res}},
#' \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.stein}},
#' \code{\link{residuals.positivestein}},
#'
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
#' model <- ipt(X, y, H, h, alpha = 0.05)
#' residuals(model)
#' @export

residuals.ipt <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' Coefficients extracted from the model object \code{improvedpreliminaryTest}
#'
#' @param object An object of class \code{improvedpreliminaryTest}.
#' @param ... Other arguments.
#'
#' @return A vector of coefficients.
#'
#' @seealso
#' \code{\link{coefficients.unres}},
#' \code{\link{coefficients.res}},
#' \code{\link{coefficients.pt}},
#' \code{\link{coefficients.st}},
#' \code{\link{coefficients.pst}},
#' \code{\link{coef.unres}},
#' \code{\link{coef.res}},
#' \code{\link{coef.pst}},
#' \code{\link{coef.st}},
#' \code{\link{coef.pst}}.
#'
#' @importFrom stats coefficients
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
#' model <- ipt(X, y, H, h, alpha = 0.05)
#' coefficients(model)
#' @export

coefficients.ipt <- function(object, ...) {
  return(object$coef)
}

#' @rdname coefficients.ipt
#' @importFrom stats coef
#' @examples
#' coef(model)
#' @export

coef.ipt <- function(object, ...) {
  return(object$coef)
}
