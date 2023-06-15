#' The positive rule Stein estimator
#'
#' This function calculates the positive-rule Stein estimatorm that is a improved version of Stein Estimator by only considering the positive part of shrinking factor
#'
#' @param X Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y Univariate quantitative response variable with dimension n_obs.
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' positiveStein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' positiveStein(X, y, H, h)
#' @export
positiveStein <- function(X, y, H, h) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)
  m <- n - p
  d <- ((q - 2) * m) / (q * (m + 2))
  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  beta <- r_est + as.numeric(1 - d / test_stat) * as.integer(test_stat > d) * (u_est - r_est)
  residuals <- y - X %*% beta
  fit <- structure(list(coefficients = beta, residuals = residuals), class = c("positiveStein"))
  fit
}



#' Predict method for Model Fits
#'
#' Predicted values based on model object.
#'
#' @param object An object of class "\code{positiveStein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{Stein}" or "\code{unrestricted}".
#' @param newdata An optional data frame in which to look for variables with which to predict.
#'  If omitted, the fitted values are used.
#' @param ... Other.
#' @seealso \code{\link{predict.positiveStein}}, \code{\link{predict.preliminaryTest}},
#' \code{\link{predict.restricted}}, \code{\link{predict.Stein}},
#' \code{\link{predict.unrestricted}}
#' \code{\link{fitted.positiveStein}}, \code{\link{fitted.preliminaryTest}},
#' \code{\link{fitted.restricted}}, \code{\link{fitted.Stein}},
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
#' model <- positiveStein(X, y, H, h)
#' fitted(model, X)
#' @export
fitted.positiveStein <- function(object, newdata, ...) {
  return(newdata %*% object$coefficients)
}

#' @rdname fitted.positiveStein
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
#' model <- positiveStein(X, y, H, h)
#' predict(model, X)
#' @export
predict.positiveStein <- function(object, newdata, ...) {
  return(newdata %*% object$coefficients)
}

#' residuals method for Model Fits
#'
#' residuals values based on model object.
#'
#' @param object An object of class "\code{positiveStein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{Stein}" or "\code{unrestricted}".
#' @param ... Other.
#' @seealso \code{\link{residuals.positiveStein}}, \code{\link{residuals.preliminaryTest}},
#' \code{\link{residuals.restricted}}, \code{\link{residuals.Stein}},
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
#' model <- positiveStein(X, y, H, h)
#' residuals(model)
#' @export

residuals.positiveStein <- function(object, ...) {
  return(object$residuals)
}

#' Extract Model Coefficients
#'
#' coef is a generic function which extracts model
#' coefficients from objects returned by modeling functions.coefficients is an alias for it.
#'
#' @param object An object of class "\code{positiveStein}", "\code{preliminaryTest}",
#' "\code{restricted}", "\code{Stein}" or "\code{unrestricted}".
#' @param ... Other.
#' @seealso \code{\link{coefficients.positiveStein}}, \code{\link{coefficients.preliminaryTest}},
#' \code{\link{coefficients.restricted}}, \code{\link{coefficients.Stein}},
#' \code{\link{coefficients.unrestricted}}
#' \code{\link{coef.positiveStein}}, \code{\link{coef.preliminaryTest}},
#' \code{\link{coef.restricted}}, \code{\link{coef.Stein}},
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
#' model <- positiveStein(X, y, H, h)
#' coefficients(model)
#' @export

coefficients.positiveStein <- function(object, ...) {
  return(object$coefficients)
}

#' @rdname coefficients.positiveStein
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
#' model <- positiveStein(X, y, H, h)
#' coef(model)
#' @export

coef.positiveStein <- function(object) {
  return(object$coefficients)
}
