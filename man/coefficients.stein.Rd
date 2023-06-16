% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/positivestein.R, R/preliminarytest.R,
%   R/restricted.R, R/stein.R, R/unrestricted.R
\name{coefficients.positivestein}
\alias{coefficients.positivestein}
\alias{coef.positivestein}
\alias{coefficients.preliminaryTest}
\alias{coef.preliminaryTest}
\alias{coefficients.restricted}
\alias{coef.restricted}
\alias{coefficients.stein}
\alias{coef.stein}
\alias{coefficients.unrestricted}
\alias{coef.unrestricted}
\title{Extract Model Coefficients}
\usage{
\method{coefficients}{positivestein}(object, ...)

\method{coef}{positivestein}(object, ...)

\method{coefficients}{preliminaryTest}(object, ...)

\method{coef}{preliminaryTest}(object, ...)

\method{coefficients}{restricted}(object, ...)

\method{coef}{restricted}(object, ...)

\method{coefficients}{stein}(object, ...)

\method{coef}{stein}(object, ...)

\method{coefficients}{unrestricted}(object, ...)

\method{coef}{unrestricted}(object, ...)
}
\arguments{
\item{object}{An object of class "\code{positivestein}", "\code{preliminaryTest}",
"\code{restricted}", "\code{stein}" or "\code{unrestricted}".}

\item{...}{Other.}
}
\description{
coef is a generic function which extracts model
coefficients from objects returned by modeling functions.coefficients is an alias for it.
}
\examples{
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- positivestein(X, y, H, h)
coefficients(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- positivestein(X, y, H, h)
coef(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- preliminaryTest(X, y, H, h, alpha = 0.05)
coefficients(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- preliminaryTest(X, y, H, h, alpha = 0.05)
coef(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- restricted(X, y, H, h)
coefficients(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- restricted(X, y, H, h)
coef(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nr = 3, nc = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- stein(X, y, H, h)
coefficients(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
p <- ncol(X)
# H beta = h
H <- matrix(c(1, 1, -1, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0), nrow = 3, ncol = p, byrow = TRUE)
h <- rep(0, nrow(H))
model <- stein(X, y, H, h)
coef(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
model <- unrestricted(X, y)
coefficients(model)
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
model <- unrestricted(X, y)
coef(model)
}
\seealso{
\code{\link{coefficients.positivestein}}, \code{\link{coefficients.preliminaryTest}},
\code{\link{coefficients.restricted}}, \code{\link{coefficients.stein}},
\code{\link{coefficients.unrestricted}}
\code{\link{coef.positivestein}}, \code{\link{coef.preliminaryTest}},
\code{\link{coef.restricted}}, \code{\link{coef.stein}},
\code{\link{coef.unrestricted}}.
}