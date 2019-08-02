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
#' @export
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
#' @param Tx A vector of case-management rates
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param X A pair of prevalence rates
fn2 <- function(params, Tx, PAR, Bfn, X) {
    A0 = params[1]
    A1 = params[2]
    D = makeD(PAR)
    PAR$A = A0
    PAR$rho = Tx[1]
    Y0 = findYeq(PAR, Bfn)
    PAR$A = A1
    PAR$rho = Tx[1]
    B1 = Bfn(PAR)
    Y1 = B1 %*% Y0
    PAR$rho = Tx[2]
    B2 = Bfn(PAR)
    Y2 = B2 %*% Y1
    return(as.numeric(sum(abs(X[1] - D %*% Y1)) + sum(abs(X[2] - D %*% Y2))))
}

#' Generic forward simulation function that can accomdate models with and
#' without exposure classes
#'
#' @param X A vector of prevalence rates
#' @param Tx A vector of case-management rates
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param cpp A toggle for using the C++ version of the forward simulation
#' algorithm
#' @export
simAR <- function(X, Tx, PAR, Bfn, cpp = T) {
    V = makeV(PAR, Bfn)
    W = makeW(PAR, Bfn)
    D = makeD(PAR)
    # outY <- matrix(nrow = length(D), ncol = length(X))
    if (PAR$In > 1) {
        ## Exposure class
        params = stats::optim(fn = fn2, par = c(0, 0.1), Tx = Tx, PAR = PAR, Bfn = Bfn, X = X)$par
        if(cpp) {
            outA <- simA2(PAR$Q, D, X, Tx, params)
        } else {
            # Set up equilibrium start values
            PAR$A = c(params[1])
            PAR$rho = Tx[1]
            Y0 = findYeq(PAR, Bfn)
            outA = PAR$A = c(params[2])
            B = Bfn(PAR)
            Y = B %*% Y0
            # outY[, 1] = Y

            # Simulate forward
            for (i in 3:length(X)) {
                PAR$rho = Tx[i - 1]
                V1 = makeV(PAR, Bfn)
                W = makeW(PAR, Bfn)
                PAR$rho = Tx[i]
                V2 = makeV(PAR, Bfn)
                outA[i - 1] = PAR$A = as.numeric((X[i] - D %*% V2 %*% V1 %*% Y)/(D %*% V2 %*% W %*% Y))
                PAR$rho = Tx[i - 1]
                B = Bfn(PAR)
                Y = B %*% Y
                # outY[, i - 1] = Y
            }
        }
    } else {
        ## No exposure class
        # Set up equilibrium start values
        outA = PAR$A = c(optimA(X[1], PAR, Bfn))
        PAR$rho = Tx[1]
        if(cpp) {
            outA = simA1(PAR$Q, D, X, Tx, outA)
        } else {
            Y = findYeq(PAR, Bfn)
            # outY[, 1] = Y
            # Simulate forward
            for (i in 2:length(X)) {
                PAR$rho = Tx[i]
                V = makeV(PAR, Bfn)
                W = makeW(PAR, Bfn)
                outA[i] = PAR$A = as.numeric((X[i] - D %*% V %*% Y)/(D %*% W %*% Y))
                PAR$rho = Tx[i]
                B = Bfn(PAR)
                Y = B %*% Y
                # outY[, i] = Y
            }
        }

    }
    # return(list(A = outA, Y = outY))
    return(outA)
}
