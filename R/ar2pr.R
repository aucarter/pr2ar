## Calculate equilibrium Y for a given A and model B
findYeq <- function(PAR, Bfn) {
    B = Bfn(PAR)
    e = eigen(B)
    first = Re(e$vectors[, 1])
    Yeq = as.vector(first/sum(first))
    return(Yeq)
}

## Calculate PR and C from input AR
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
