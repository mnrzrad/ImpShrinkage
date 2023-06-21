#' Test-Statistics
#'
#' This function calculates the test statistics, assuming
#' \eqn{\mathcal{H}_0: H \beta = h}. When the error has a normal distribution,
#' it is defined as
#' \deqn{\mathcal{L} = \frac{(H\hat{\beta}^{U}-h)^{\top}(H(X^{\top}X)^{-1}
#' H^{\top})^{-1}(H\hat{\beta}^{U}-h) }{q s^2_{unr}}}
#' and when the error has a non-normal distribution, as
#' \deqn{\mathcal{L} = \frac{(H\hat{\beta}^{U}-h)^{\top}(H(X^{\top}X)^{-1}
#' H^{\top})^{-1}(H\hat{\beta}^{U}-h) }{s^2_{unr}}}
#' where
#' \itemize{
#'   \item \eqn{\hat{\beta}^{U}} is the unrestricted estimator; See \code{\link{unrReg}}.
#'   \item \eqn{q} is the number of restrictions, i.e., the number of rows of known matrix
#'   \eqn{H};
#'   \item \eqn{s^2_{unr}} is the corresponding unrestricted estimator of
#'   \eqn{\sigma^2}.
#' }
#'
#'
#' @param X Matrix with input observations, of dimension \code{n} x \code{p};
#' each row is an observation vector.
#' @param y Vector with response observations of size \code{n}.
#' @param H A given \code{q} x \code{p} matrix.
#' @param h A given \code{q} x \code{1} vector.
#' @param is_error_normal logical value indicating whether the errors follow a
#' normal distribution. If \code{is_error_normal} is \code{TRUE}, the distribution
#' of the test statistics for the null hypothesis is the F distribution
#' (\code{\link[stats]{FDist}}).On the other hand, if the errors have a
#' non-normal distribution, the asymptotic distribution of the test statistics
#' is the \eqn{\chi^2} distribution (\code{\link[stats]{Chisquare}}). By default,
#' \code{is_error_normal} is set to \code{FALSE}.
#'
#' @return The value of the test statistic.
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
#' simulated_data <- simdata(n = n_obs, p_vars, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' p <- ncol(X)
#' # H beta = h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' teststat(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, nrow(H))
#' teststat(X, y, H, h)
#'
#' data(cement)
#' X <- as.matrix(cbind(1, cement[, 1:4]))
#' y <- cement$y
#' # Based on Kaciranlar et al. (1999)
#' H <- matrix(c(0, 1, -1, 1, 0), nrow = 1, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' teststat(X, y, H, h)
#' # Based on Kibria (2005)
#' H <- matrix(c(0, 1, -1, 1, 0, 0, 0, 1, -1, -1, 0, 1, -1, 0, -1), nrow = 3, ncol = 5, byrow = TRUE)
#' h <- rep(0, nrow(H))
#' teststat(X, y, H, h)
#' @export
teststat <- function(X, y, H, h, is_error_normal = FALSE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- n - p
  q <- nrow(H)
  u_est <- unrReg(X, y) # unrestricted estimator
  s2 <- (t(y - X %*% u_est$coef) %*% (y - X %*% u_est$coef)) / m
  C <- t(X) %*% X
  diff <- (H %*% u_est$coef - h)
  if (!is_error_normal) {
    as.numeric((t(diff) %*% solve(H %*% solve(C) %*% t(H)) %*% diff) / (q * s2))
  } else {
    as.numeric((t(diff) %*% solve(H %*% solve(C) %*% t(H)) %*% diff) / (s2))
  }
}
