#' The restricted estimator
#'
#' This function calculates the restricted estimator using the formula:
#' \deqn{\hat{\beta}^{R} = \hat{\beta}^{U} - (X^{\top}X)^{-1}H^{\top}
#' (H(X^{\top}X)^{-1}H^{\top})^{-1}(H\hat{\beta}^{U}-h)}
#' where \eqn{H\beta = h} represents a subspace of the parameter space induced
#' by the non-sample information. Here, \eqn{H} is a known \eqn{q \times p}
#' matrix, and \eqn{h} is a known \eqn{q}-vector.
#'
#'
#' The corresponding unrestricted estimator of \eqn{\sigma^2} is
#' \deqn{s^2 = \frac{1}{n-p}(y-X\hat{\beta}^{R})^{\top}(y - X\hat{\beta}^{R})}
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#'
#' @return A vector of regression coefficients
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0),
#'   nrow = 3,
#'   ncol = p, byrow = TRUE
#' )
#' h <- rep(0, nrow(H))
#' restricted(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0),
#'   nrow = 3,
#'   ncol = p, byrow = TRUE
#' )
#' h <- rep(1, nrow(H))
#' restricted(X, y, H, h)
#'
#'
#' data(cement)
#' X <- as.matrix(cbind(1,cement[,1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0,1,-1,1,0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' restricted(X, y, H, h)
#'
#' H <- matrix(c(0,1,-1,1,0,0, 0, 1, -1, -1,0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' restricted(X, y, H, h)
#' @export

restricted <- function(X, y, H, h) {
  d <- dim(X)
  u_est <- unrestricted(X, y) # unrestricted estimator
  C <- t(X) %*% X
  beta <- u_est$coef - solve(C) %*% t(H) %*% solve(H %*% solve(C) %*% t(H)) %*% (H %*% u_est$coef - h)
  residuals <- (y - X %*% beta)[, 1]
  n <- dim(X)[1]
  p <- dim(X)[2]
  s2 <- sum(residuals^2) / (n - p)
  fittedValues <- (X %*% beta)[, 1]
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2, fitted.value = fittedValues), class = c("restricted"))
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
#' model <- restricted(X, y, H, h)
#' fitted(model)
#' @export
fitted.restricted <- function(object, ...) {
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
#' model <- restricted(X, y, H, h)
#' predict(model, X)
#' @export
predict.restricted <- function(object, newdata, ...) {
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
#' model <- restricted(X, y, H, h)
#' residuals(model)
#' @export

residuals.restricted <- function(object, ...) {
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
#' model <- restricted(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.restricted <- function(object, ...) {
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
#' model <- restricted(X, y, H, h)
#' coef(model)
#' @export

coef.restricted <- function(object, ...) {
  return(object$coef)
}
