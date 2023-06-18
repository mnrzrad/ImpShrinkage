#' The improved preliminary test estimator
#'
#' This function calculates the improved preliminary test estimator that is
#' calculated by
#' If the error has a normal distribution:#'
#' \deqn{\hat{\beta}^{ipt}= \hat{\beta}^{PT} - d (\hat{\beta}^{U} - \hat{\beta}^{R}) \mathcal{L}^{-1} I(\mathcal{L} > F_{q,m}(\alpha)) }
# If the error has a non-normal distribution:
#' \deqn{\hat{\beta}^{ipt}= \hat{\beta}^{PT} - d (\hat{\beta}^{U} - \hat{\beta}^{R}) \mathcal{L}^{-1} I(\mathcal{L} > \chi^2_{q}(\alpha)) }
#``
#' where \eqn{I(A)} denotes an indicator function and
#' \itemize{
#'   \item \eqn{\hat{\beta}^{PT}}: the \code{\link{preliminaryTest}} estimator
#'   \item \eqn{\hat{\beta}^{U}}: the \code{\link{unrestricted}} estimator
#'   \item \eqn{\hat{\beta}^{R}}: the \code{\link{restricted}} estimator
#'   \item \eqn{\mathcal{L}}: the \code{\link{test_staistics}}
#'   \item \eqn{F_{q,m}(\alpha)}: the upper \eqn{\alpha} level critical value of \eqn{F}-distribution with \eqn{(q,n-p)} degrees of freedom, calculated using \code{\link[stats]{qf}}.
#'   \item \eqn{\chi^2_{q}(\alpha)}: the upper \eqn{\alpha} level critical value of \eqn{\chi^2}-distribution with \eqn{q} degree of freedom, calculated using \code{\link[stats]{qchisq}}.
#'   \item \eqn{d}: the shrinkage factor
#'   \item \eqn{\alpha}: the significance level.
#' }
#'
#'
#' #'The corresponding unrestricted estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{iPT})^{\top}(y - X\hat{\beta}^{iPT})}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#' @param alpha  A given significance level
#' @param d An optional parameter. If not provided (or set to \code{NULL}), it will be
#' calculated using \eqn{\frac{{(q - 2) \cdot (n - p}}{{q \cdot (n - p + 2)}}}
#' @param normal_error logical value indicating whether the errors follow a
#' normal distribution. #'If \code{normal_error} is \code{TRUE}, the distribution
#' of the test statistics for the null hypothesis is F distribution,
#' \code{\link[stats]{FDist}}. On the other hand, if the errors have a
#' non-normal distribution, the asymptotic distribution of the test statistics
#' is \eqn{\chi^2} distribution, \code{\link[stats]{Chisquare}}. By default,
#' \code{normal_error} is set to \code{FALSE}.
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
#' improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' # H beta != h
#' p <- ncol(X)
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#'
#' @importFrom stats qf
#' @export

improvedpreliminaryTest <- function(X, y, H, h, alpha, d = NULL, normal_error = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  q <- nrow(H)
  m <- n - p
  if (is.null(d)) {
    d <- ((q - 2) * m) / (q * (m + 2))
  }
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  pt_est <- preliminaryTest(X, y, H, h, alpha)
  test_stat <- test_statistics(X, y, H, h, normal_error = normal_error)
  if (!normal_error) {
    threshold <- stats::qf(1 - alpha, q, n - p)
  } else {
    threshold <- stats::qchisq(1 - alpha, q)
  }
  threshold <- stats::qf(1 - alpha, q, n - p)
  beta <- pt_est$coef - d * ((u_est$coef - r_est$coef) / test_stat) * as.integer(test_stat < threshold)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2, fitted.value = fittedValues), class = c("improvedpreliminaryTest"))
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#' fitted(model, X)
#' @export
fitted.improvedpreliminaryTest <- function(object, newdata, ...) {
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
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#' predict(model, X)
#' @export
predict.improvedpreliminaryTest <- function(object, newdata, ...) {
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
#' model <- improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#' residuals(model)
#' @export

residuals.improvedpreliminaryTest <- function(object, ...) {
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
#' model <- improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#' coefficients(model)
#' @export

coefficients.improvedpreliminaryTest <- function(object, ...) {
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
#' model <- improvedpreliminaryTest(X, y, H, h, alpha = 0.05)
#' coef(model)
#' @export

coef.improvedpreliminaryTest <- function(object, ...) {
  return(object$coef)
}
