## Generic make V and make W function for any input B function
makeV = function(PAR, Bfn) {
    PAR$A = 0
    Bfn(PAR)
}

makeW = function(PAR, Bfn) {
    V = makeV(PAR, Bfn)
    (Bfn(PAR) - V)/PAR$A
}

## Make D to map Y to X
makeD <- function(PAR) {
    In = PAR$In
    Cn = PAR$Cn
    D = c(0)
    if (In > 0) {
        D = c(D, 0, rep(1, In * 2))
    } else {
        D = c(D, 1)
    }
    D = c(D, rep(0, Cn))
    return(D)
}

## Simple case
makeB <- function(PAR) {
    A = PAR$A
    Q = PAR$Q
    cbind(c(1 - A, A), c((1 - A) * (1 - Q), A * (1 - Q) + Q))
    
}

## Drugs
makeBdrugs <- function(PAR) {
    A = PAR$A
    Q = PAR$Q
    d = PAR$d
    rho = PAR$rho
    cbind(c((1 - d) * (1 - A), (1 - d) * (1 - rho) * A, d + (1 - d) * rho * A), c((1 - d) * (1 - A * 
        rho) * (1 - Q), (1 - d) * (1 - A * rho) * Q, d + (1 - d) * rho * A), c(1 - d, 0, d))
}

## Age
makeBage <- function(PAR) {
    A = PAR$A
    Q = PAR$Q
    In = PAR$In
    S = c(S = 1 - A, E0 = A, E = rep(0, In), I = rep(0, In))
    E0 = c(0, 0, A, rep(0, In - 1), (1 - A), rep(0, In - 1))
    M = cbind(S, E0)
    
    Ei = E0
    for (i in 1:In) M = cbind(M, Ei)
    
    for (i in 1:(In - 1)) {
        Ii = c((1 - A) * (1 - Q), A * (1 - Q), rep(0, i), A * Q, rep(0, In - i - 1), rep(0, i), (1 - 
            A) * Q, rep(0, In - i - 1))
        M = cbind(M, Ii)
    }
    
    Im = c((1 - A) * (1 - Q), A * (1 - Q), rep(0, In - 1), A * Q, rep(0, In - 1), (1 - A) * Q)
    M = cbind(M, Im)
    
    M
}

## Drugs + Age
makeBdrugs_age = function(PAR) {
    A = PAR$A
    Q = PAR$Q
    rho = PAR$rho
    In = PAR$In
    Cn = PAR$Cn
    S = c(S = 1 - A, E0 = A, E = rep(0, In), I = rep(0, In), C = rep(0, Cn))
    E0 = c(0, 0, (1 - rho) * A, rep(0, In - 1), (1 - rho) * (1 - A), rep(0, In - 1), rho, rep(0, Cn - 
        1))
    M = cbind(S, E0)
    
    Ei = E0
    for (i in 1:In) M = cbind(M, Ei)
    
    for (i in 1:(In - 1)) {
        Ii = c((1 - A) * (1 - Q), A * (1 - Q), rep(0, i), A * Q, rep(0, In - i - 1), rep(0, i), (1 - 
            A) * Q, rep(0, In - i - 1), rep(0, Cn))
        M = cbind(M, Ii)
    }
    
    Im = c((1 - A) * (1 - Q), A * (1 - Q), rep(0, In - 1), A * Q, rep(0, In - 1), (1 - A) * Q, rep(0, 
        Cn))
    M = cbind(M, Im)
    
    if (Cn == 2) {
        C1 = c(0, rep(0, 2 * In + 1), 0, 1)
        M = cbind(M, C1)
        C2 = c(1, rep(0, 2 * In + 1), 0, 0)
        M = cbind(M, C2)
    }
    M
}
