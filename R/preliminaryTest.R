preliminaryTest <- function(X,y, H, h, q, alpha){
    n <- dim(X)[1]
    p <- dim(X)[2]
    u_est <- unrestricted(X, y)
    r_est <- restricted(X, y, H, h)
    test_stat <- test_statistics(X, y, H, h, q)
    threshold <- qf(1-alpha,q, n-p)
    u_est - (u_est-r_est) * integer(test_stat < threshold)
}


