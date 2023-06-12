generate_H_and_h <- function(p) {
    H <- matrix(0, nrow = p - 1, ncol = p)
    for (i in 1:(p - 1)) {
        H[i, i] <- 1
        H[i, i+1] <- -1
    }
    h <- rep(0, p-1)
    return(list(H = H, h = h))
}

