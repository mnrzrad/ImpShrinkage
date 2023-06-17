#' The improved preliminary test estimator
#'
#' This function calculates the improved preliminary test for the null hypothesis: H.beta = h
#'
#' @param X  Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y  Univariate quantitative response variable with dimension n_obs
#' @param H  A given q_restr x p_vars matrix.
#' @param h  A given q_restr x 1 vector.
#' @param alpha  A given significance level
#' @param d An optional parameter. If not provided (or set to NULL), it will be calculated using \eqn{\frac{{(q_restr - 2) \cdot (n_obs - p_vars}}{{q_restr \cdot (n_obs - p_vars + 2)}}}
#' .
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

improvedpreliminaryTest <- function(X, y, H, h, alpha, d = NULL) {
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
  test_stat <- test_statistics(X, y, H, h)
  threshold <- stats::qf(1 - alpha, q, n - p)
  beta <- pt_est$coef - d * ((u_est$coef - r_est$coef) / test_stat) * as.integer(test_stat < threshold)
  residuals <- (y - X %*% beta)[, 1]
  s2 <- sum(residuals^2) / (n - p)
  fit <- structure(list(coef = beta, residuals = residuals, s2 = s2), class = c("improvedpreliminaryTest"))
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
#' model <- improvecpreliminaryTest(X, y, H, h, alpha = 0.05)
#' fitted(model, X)
#' @export
fitted.preliminaryTest <- function(object, newdata, ...) {
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
#' model <- imrpovedpreliminaryTest(X, y, H, h, alpha = 0.05)
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
