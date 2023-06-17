#' The Stein estimator
#'
#' This function can be used to calculate the Stein estimator
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs
#' @param H A given q_restr x p_vars matrix.
#' @param h A given q_restr x 1 vector.
#'
#' @return A vector of regression coefficients
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
#' stein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' stein(X, y, H, h)
#' @export

stein <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  beta <- u_est$coef - d * (u_est$coef - r_est$coef) / test_stat
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2), class = c("stein"))
  fit
}


#' Predict method for Model Fits
#'
#' Predicted values based on model object.
#'
#' @param object An object of class "\code{stein}, "\code{preliminaryTest}",
#' "\code{restricted}", "\code{positivestein}"" or "\code{unrestricted}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other.
#' @seealso \code{\link{predict.positivestein}}, \code{\link{predict.preliminaryTest}},
#' \code{\link{predict.restricted}}, \code{\link{predict.stein}},
#' \code{\link{predict.unrestricted}}
#' \code{\link{fitted.positivestein}}, \code{\link{fitted.preliminaryTest}},
#' \code{\link{fitted.restricted}}, \code{\link{fitted.stein}},
#' \code{\link{fitted.unrestricted}}.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- stein(X, y, H, h)
#' fitted(model, X)
#' @export
fitted.stein <- function(object, newdata, ...) {
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- stein(X, y, H, h)
#' predict(model, X)
#' @export
predict.stein <- function(object, newdata, ...) {
  return((newdata %*% object$coef)[, 1])
}

#' residuals method for Model Fits
#'
#' residuals values based on model object.
#'
#' @param object An object of class "\code{positivestein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{stein}" or "\code{unrestricted}".
#' @param ... Other.
#' @seealso \code{\link{residuals.positivestein}}, \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.restricted}}, \code{\link{residuals.stein}},
#' \code{\link{residuals.unrestricted}}.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- stein(X, y, H, h)
#' residuals(model)
#' @export

residuals.stein <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' coef is a generic function which extracts model
#' coefficients from objects returned by modeling functions.coefficients is an alias for it.
#'
#' @param object An object of class "\code{positivestein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{stein}" or "\code{unrestricted}".
#' @param ... Other.
#' @seealso \code{\link{coefficients.positivestein}}, \code{\link{coefficients.preliminaryTest}},
#' \code{\link{coefficients.restricted}}, \code{\link{coefficients.stein}},
#' \code{\link{coefficients.unrestricted}}
#' \code{\link{coef.positivestein}}, \code{\link{coef.preliminaryTest}},
#' \code{\link{coef.restricted}}, \code{\link{coef.stein}},
#' \code{\link{coef.unrestricted}}.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- stein(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.stein <- function(object, ...) {
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
#' model <- stein(X, y, H, h)
#' coef(model)
#' @export

coef.stein <- function(object, ...) {
  return(object$coef)
}
