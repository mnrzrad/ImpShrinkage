#' The positive rule Stein estimator
#'
#' positiveStein
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
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(0, q)
#' positiveStein(X, y, H, h)
#'
#' # H beta != h
#' H <- matrix(c(1,1,-1,0,0,1,0,1,0,-1,0,0,0,1,0), nr = 3, nc = p, byrow = TRUE)
#' h <- rep(1, q)
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
  test_stat <- test_statistics(X, y, H, h, q)
  return(r_est + as.numeric(1 - d / test_stat) * as.integer(test_stat > d) * (u_est - r_est))
}


