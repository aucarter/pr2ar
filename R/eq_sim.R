
## Pair of function for finding the maximum value of X and associated A
AR2PRopt <- function(A, PAR, Bfn) {
    AR2PR(A, PAR, Bfn)$X
}

findAXmax <- function(PAR, Bfn) {
    opt <- stats::optimize(AR2PRopt, interval = c(0, 1), PAR, Bfn, maximum = T)
    A = opt$maximum
    X = opt$objective
    return(list(A = A, X = X))
}

## Calculate equilibrium AR for a PR time-series
PR2AReq <- function(X, PAR, Bfn, showMessages = F) {
    A1X <- AR2PR(1, PAR, Bfn)
    max <- findAXmax(PAR, Bfn)
    mono <- abs(A1X$X - max$X) < 0.001
    if (showMessages) {
        print(paste0("max A: ", round(max$A, 3), "; max X: ", round(max$X, 3)))
    }
    outA <- c()
    outY <- matrix(nrow = nrow(Bfn(PAR)), ncol = length(X))
    if (!mono) {
        outA2 <- c()
        outY2 <- matrix(nrow = nrow(Bfn(PAR)), ncol = length(X))
    }
    for (Xi in X) {
        if (Xi > max$X) {
            outA <- c(outA, NA)
            outA2 <- c(outA2, NA)
            next
        }
        A = stats::optimize(fn, c(0, max$A), PAR = PAR, X = Xi, Bfn = Bfn)$minimum
        outA = c(outA, A)
        PAR$A = A
        Y = findYeq(PAR, Bfn)
        outY[, which(X == Xi)] = Y
        if (!mono) {
            if (Xi > A1X$X) {
                A = stats::optimize(fn, c(max$A, 1), PAR = PAR, X = Xi, Bfn = Bfn)$minimum
                outA2 <- c(outA2, A)
                PAR$A = A
                Y = findYeq(PAR, Bfn)
                outY2[, which(X == Xi)] <- Y
            } else {
                outA2 <- c(outA2, NA)
            }
        }
    }
    if (mono) {
        outList = list(A = outA, Y = outY)
    } else {
        outList = list(A = outA, Y = outY, A2 = outA2, Y2 = outY2)
    }
    return(outList)
}
