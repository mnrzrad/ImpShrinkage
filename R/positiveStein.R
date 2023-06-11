positiveStein <- function(X,y, H, h, q){
    n <- dim(X)[1]
    p <- dim(X)[2]
    m <- n - p
    d <- ((q-2)*m)/(q*(m+2))
    u_est <- unrestricted(X, y)
    r_est <- restricted(X, y, H, h)
    test_stat <- test_statistics(X, y, H, h, q)
    threshold <- qf(1-alpha,q, n-p)
    r_est + (1 - d/test_stat)* int(test_stat > d) * (u_est-r_est)
}


