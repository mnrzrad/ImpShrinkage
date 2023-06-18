#' The preliminary test estimator
#'
#' This function calculates the preliminary test using the formula:
#' If the error has a normal distribution:
#' \deqn{\hat{\beta}^{PT}=\hat{\beta}^{U} - (\hat{\beta}^{U} - \hat{\beta}^{R}) I(\mathcal{L} \le F_{q,m}(\alpha))}
#' If the error has a non-normal distribution:
#' \deqn{\hat{\beta}^{PT}=\hat{\beta}^{U} - (\hat{\beta}^{U} - \hat{\beta}^{R}) I(\mathcal{L} \le \chi^2_{q}(\alpha))}
#' #' where \eqn{I(A)} denotes an indicator function and
#'\itemize{
#'   \item \eqn{\hat{\beta}^{U}}: the \code{\link{unrestricted}} estimator.
#'   \item \eqn{\hat{\beta}^{R}}: the \code{\link{restricted}} estimator.
#'   \item \eqn{\mathcal{L}}: the \code{\link{test_statistics}}.
#'   \item \eqn{F_{q,m}(\alpha)}: the upper \eqn{\alpha} level critical value of \eqn{F}-distribution with \eqn{(q,n-p)} degrees of freedom, calculated using \code{\link[stats]{qf}}.
#'   \item \eqn{\chi^2_{q}(\alpha)}: the upper \eqn{\alpha} level critical value of \eqn{\chi^2}-distribution with \eqn{q} degree of freedom, calculated using \code{\link[stats]{qchisq}}.
#'   \item \eqn{\alpha}: the significance level.
#' }
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#' @param alpha  A given significance level.
#' @param normal_error logical value indicating whether the errors follow a normal distribution. #'If \code{normal_error} is \code{TRUE}, the distribution of the test statistics for the null hypothesis is F distribution, \code{\link[stats]{FDist}}.
#'  On the other hand, if the errors have a non-normal distribution, the asymptotic distribution of the test statistics is \eqn{\chi^2} distribution, \code{\link[stats]{Chisquare}}. By default, \code{normal_error} is set to \code{FALSE}.
#'
#' @return  A vector of regression coefficients
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' # H beta != h
#' p <- ncol(X)
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' preliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' @importFrom stats qf
#' @export

preliminaryTest <- function(X, y, H, h, alpha, normal_error = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h, normal_error = normal_error)
  if (!normal_error) {
    threshold <- stats::qf(1 - alpha, q, n - p)
  } else {
    threshold <- stats::qchisq(1 - alpha, q)
  }
  beta <- u_est$coef - (u_est$coef - r_est$coef) * as.integer(test_stat < threshold)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2, fitted.value = fittedValues), class = c("preliminaryTest"))
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' fitted(model)
#' @export
fitted.preliminaryTest <- function(object, ...) {
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' predict(model, X)
#' @export
predict.preliminaryTest <- function(object, newdata, ...) {
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' residuals(model)
#' @export

residuals.preliminaryTest <- function(object, ...) {
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' coefficients(model)
#' @export

coefficients.preliminaryTest <- function(object, ...) {
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' coef(model)
#' @export

coef.preliminaryTest <- function(object, ...) {
  return(object$coef)
}
