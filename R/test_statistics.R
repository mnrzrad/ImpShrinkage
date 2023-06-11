test_statistics <- function(X, y, H, h, q){
    n <- dim(X)[1]
    p <- dim(X)[2]
    u_est <- unrestricted(X, y) #unrestricted estimator
    s2 <- (t(y - X %*% u_est) %*% (y - X %*% u_est))/(n-p)
    C <- t(X) %*% X
    diff <- (H %*% u_est - h)
    (t(diff) %*% solve(H %*% solve(C) %*% t(H)) %*% diff) / (q * s2)
}
