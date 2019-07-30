#' Calculates the equilibrium Y for a given attack-rate and model
#'
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
findYeq <- function(PAR, Bfn) {
    B = Bfn(PAR)
    e = eigen(B)
    first = Re(e$vectors[, 1])
    Yeq = as.vector(first/sum(first))
    return(Yeq)
}

#' Calculate the prevalence rate and treatment rate for a given attack-rate,
#' model, and set of model parameters
#'
#' @param A An attack-rate
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
AR2PR <- function(A, PAR, Bfn) {
    D = makeD(PAR)
    PR = c()
    for (Ai in A) {
        PAR$A = Ai
        eq = findYeq(PAR, Bfn)
        PR = c(PR, D %*% eq)
        if (PAR$Cn > 0) {
            Ci = sum(eq[(length(eq) - 1:PAR$Cn) + 1])
            if (Ai == A[1]) {
                C = c(Ci)
            } else {
                C = c(C, Ci)
            }
        }
    }
    if (length(eq) == 3) {
        out.list = list(X = PR, C = C)
    } else {
        out.list = list(X = PR)
    }
    return(out.list)
}
