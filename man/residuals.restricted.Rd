% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resReg.R
\name{residuals.restricted}
\alias{residuals.restricted}
\title{Extract Model Residuals}
\usage{
\method{residuals}{restricted}(object, ...)
}
\arguments{
\item{object}{An object of class \code{restricted}.}

\item{...}{Other arguments.}
}
\value{
A vector of residuals.

\code{\link{residuals.unrestricted}},
\code{\link{residuals.preliminaryTest}},
\code{\link{residuals.improvedpreliminaryTest}},
\code{\link{residuals.stein}},
\code{\link{residuals.positivestein}}.
}
\description{
Residuals values based on model object \code{restricted}.
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
model <- resReg(X, y, H, h)
residuals(model)
}
