#' The preliminary test estimator
#'
#' This function calculates the improved estimator after performing the preliminary test for the null hipothesis: H.beta = h
#'
#' @param X  Matrix with input observations, of dimension n_obs x p_vars; each row is an observation vector.
#' @param y  Univariate quantitative response variable with dimension n_obs
#' @param H  A given q_restr x p_vars matrix.
#' @param h  A given q_restr x 1 vector.
#' @param alpha  A given significance level
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

preliminaryTest <- function(X, y, H, h, alpha) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- nrow(H)

  u_est <- unrestricted(X, y)
  r_est <- restricted(X, y, H, h)
  test_stat <- test_statistics(X, y, H, h)
  threshold <- stats::qf(1 - alpha, q, n - p)
  beta <- u_est - (u_est - r_est) * as.integer(test_stat < threshold)
  residuals <- y - X %*% beta
  fit <- structure(list(coef = beta, residuals = residuals), class = c("preliminaryTest"))
  fit
}




#' @rdname fitted.positiveStein
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
#' fitted(model, X)
#' @export
fitted.preliminaryTest <- function(object, newdata, ...) {
    return(newdata %*% object$coef)
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
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' predict(model, X)
#' @export
predict.preliminaryTest <- function(object, newdata, ...) {
    return(newdata %*% object$coef)
}

#' @rdname residuals.positiveStein
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
#' model <- preliminaryTest(X, y, H, h, alpha = 0.05)
#' coefficients(model)
#' @export

coefficients.preliminaryTest <- function(object, ...) {
    return(object$coef)
}

#' @rdname coefficients.positiveStein
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
