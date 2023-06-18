#' The positive-rule Stein estimator
#'
#' This function calculates the positive-rule Stein estimator that is an improved
#' version of Stein Estimator by only considering the positive part of
#' shrinking factor. It is calculated by
#' \deqn{\hat{\beta}^{S+}= \hat{\beta}^{S} + (1 + d \mathcal{L}^{-1}) I(\mathcal{L} > d) (\hat{\beta}^{U} - \hat{\beta}^{R})}
#' where \eqn{I(A)} denotes an indicator function and
#' \itemize{
#'   \item \eqn{\hat{\beta}^{S}}: the \code{\link{stein}} estimator
#'   \item \eqn{\hat{\beta}^{U}}: the \code{\link{unrestricted}} estimator
#'   \item \eqn{\hat{\beta}^{R}}: the \code{\link{restricted}} estimator
#'   \item \eqn{\mathcal{L}}: the \code{\link{test_statistics}}
#'   \item \eqn{d}: the shrinkage factor
#' }
#'
#' The corresponding unrestricted estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{S+})^{\top}(y - X\hat{\beta}^{S+})}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
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
#'
#' @return A vector of regression coefficients
#'
#' @references
#'  Saleh, A. K. Md. Ehsanes. (2006). \emph{Theory of Preliminary Test
#'   and Stein‐Type Estimation With Applications}, Wiley.
#'
#'  ref
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
#' @export
positivestein <- function(X, y, H, h, d = NULL, normal_error = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h, normal_error = normal_error)
  beta <- r_est$coef + as.numeric(1 - d / test_stat) * as.integer(test_stat > d) * (u_est$coef - r_est$coef)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2, fitted.value = fittedValues), class = c("positivestein"))
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
#' model <- positivestein(X, y, H, h)
#' fitted(model)
#' @export
fitted.positivestein <- function(object, ...) {
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
#' model <- positivestein(X, y, H, h)
#' predict(model, X)
#' @export
predict.positivestein <- function(object, newdata, ...) {
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
#' model <- positivestein(X, y, H, h)
#' residuals(model)
#' @export

residuals.positivestein <- function(object, ...) {
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
#' model <- positivestein(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.positivestein <- function(object, ...) {
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
#' model <- positivestein(X, y, H, h)
#' coef(model)
#' @export

coef.positivestein <- function(object, ...) {
  return(object$coef)
}
