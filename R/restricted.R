restrcited <- function(X,y, H, h){
    u_est <- unrestricted(X, y) #unrestricted estimator
    C <- t(X) %*% X
    u_est - solve(C) %*% t(H) %*% solve(H %*% solve(C) %*% t(H)) %*% (H %*% u_est - h)
}
