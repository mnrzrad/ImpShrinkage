simulate_data <- function(n, p, beta){
    X <- matrix(rnorm(n*p), nr = n, nc = p)
    y <- X %*% beta + rnorm(n)
}
