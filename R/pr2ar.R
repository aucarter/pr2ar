## Set up PR2AR that can take any model and return equilibrium or forward simulation values
PR2AR <- function(X, PAR, eq = F, showMessages = F, onlyAR = F) {
    if (any(is.na(X))) {
        if (onlyAR) {
            return(rep(NA, length(X)))
        } else {
            return(list(A = rep(NA, length(X)), Y = matrix(NA, nrow = length(makeD(PAR)), ncol = length(X))))
        }
    }
    if (PAR$In > 1 & PAR$Cn > 0) {
        Bfn = makeBdrugs_age
    } else if (PAR$In > 1) {
        Bfn = makeBage
    } else if (PAR$Cn > 0) {
        Bfn = makeBdrugs
    } else {
        Bfn = makeB
    }
    if (eq) {
        outList <- PR2AReq(X, PAR, Bfn, showMessages)
    } else {
        outList <- simAR(X, PAR, Bfn)
    }
    if (onlyAR) {
        return(outList$A)
    } else {
        return(outList)
    }
}
