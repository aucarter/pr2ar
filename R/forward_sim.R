#' Helper objective function for finding the optimal attack-rate for a given
#' prevalence rate
#'
#' @param A An attack-rate
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param X A prevalence rate
fn <- function(A, PAR, Bfn, X) {
    PR = AR2PR(A, PAR, Bfn)$X
    return(abs(X - PR))
}
#' A function for finding the optimal attack-rate for a given prevalence rate
#'
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param X A prevalence rate
## Find optimal A
optimA <- function(X, PAR, Bfn) {
    A <- c()
    for(Xi in X) {
        Ai <- stats::optimize(fn, c(0, 1), PAR = PAR, X = Xi, Bfn = Bfn)$minimum
        A <- c(A, Ai)
    }
    return(A)
}

#' Helper objective function for finding a pair of optimal attack-rates for a
#' given pair of starting prevalence rates when the model has an exposure class
#'
#' @param params A pair of attack-rates
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param X A pair of prevalence rates
fn2 <- function(params, PAR, Bfn, X) {
    A1 = params[1]
    A2 = params[2]
    D = makeD(PAR)
    PAR$A = A1
    Y0 = findYeq(PAR, Bfn)
    PAR$A = A2
    B2 = Bfn(PAR)
    Y1 = B2 %*% Y0
    Y2 = B2 %*% Y1
    return(as.numeric(sum(abs(X[1] - D %*% Y1)) + sum(abs(X[2] - D %*% Y2))))
}

#' Generic forward simulation function that can accomdate models with and
#' without exposure classes
#'
#' @param X A vector of prevalence rates
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
simAR <- function(X, PAR, Bfn) {
    V = makeV(PAR, Bfn)
    W = makeW(PAR, Bfn)
    D = makeD(PAR)
    outY <- matrix(nrow = length(D), ncol = length(X))
    if (PAR$In > 1) {
        params = stats::optim(fn = fn2, par = c(0, 0.1), PAR = PAR, Bfn = Bfn, X = X)$par
        # Set up equilibrium start values
        outA = PAR$A = c(params[1])
        Y0 = findYeq(PAR, Bfn)
        outA[2] = PAR$A = c(params[2])
        B = Bfn(PAR)
        Y = B %*% Y0
        outY[, 1] = Y

        # Simulate forward
        for (i in 3:length(X)) {
            outA[i] = PAR$A = as.numeric((X[i] - D %*% V %*% V %*% Y)/(D %*% V %*% W %*% Y))
            B = Bfn(PAR)
            Y = B %*% Y
            outY[, i - 1] = Y
        }
    } else {
        # Set up equilibrium start values
        forwardA = PAR$A = optimA(X[1], PAR, Bfn)
        Y = findYeq(PAR, Bfn)
        outY[, 1] = Y
        # Simulate forward
        for (i in 2:length(X)) {
            forwardA[i] = PAR$A = as.numeric((X[i] - D %*% V %*% Y)/(D %*% W %*% Y))
            B = Bfn(PAR)
            Y = B %*% Y
            outY[, i] = Y
        }
    }
    return(list(A = outA, Y = outY))
}
