unrestrcited <- function(X,y){
    solve(t(X)%*%X) %*% t(X) %*% y
}
