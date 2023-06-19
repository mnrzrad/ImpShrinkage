% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unrestricted.R
\name{residuals.unrestricted}
\alias{residuals.unrestricted}
\title{residuals method for Model Fits}
\usage{
\method{residuals}{unrestricted}(object, ...)
}
\arguments{
\item{object}{An object of class "\code{unrestricted}".}

\item{...}{Other.}
}
\description{
residuals values based on model object.
}
\examples{
n_obs <- 100
p_vars <- 5
beta <- c(2, 1, 3, 0, 5)
simulated_data <- simdata(n = n_obs, p = p_vars, beta)
X <- simulated_data$X
y <- simulated_data$y
model <- unrestricted(X, y)
residuals(model)
}
\seealso{
\code{\link{residuals.positivestein}}, \code{\link{residuals.preliminaryTest}},
\code{\link{residuals.restricted}}, \code{\link{residuals.stein}},
\code{\link{residuals.unrestricted}}, \code{\link{residuals.improvedpreliminaryTest}}.
}