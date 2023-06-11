#' Test-Statistics
#'
#' Stein
#'
#' @param x scaler.
#' @param q scaler.
#' @param n number.
#' @param G scaler.
#'
#' @references ref
#' @examples
#' n <- 100
#' p <- 5
#' beta <- c(2, 1, 3, 0, 5)
#' simulated_data <- simulate_data(n, p, beta)
#' X <- simulated_data$X
#' y <- simulated_data$y
#' # H beta = h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(0, q)
#' test_statistics(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nrow = 3, ncol = p, byrow = TRUE)
#' h <- rep(1, q)
#' test_statistics(X, y, H, h)
#' @export




test_statistics <- function(X, y, H, h){
    n <- dim(X)[1]
    p <- dim(X)[2]
    m <- n - p
    q <- nrow(H)
    u_est <- unrestricted(X, y) #unrestricted estimator
    s2 <- (t(y - X %*% u_est) %*% (y - X %*% u_est))/m
    C <- t(X) %*% X
    diff <- (H %*% u_est - h)
    as.numeric((t(diff) %*% solve(H %*% solve(C) %*% t(H)) %*% diff) / (q * s2))
}
