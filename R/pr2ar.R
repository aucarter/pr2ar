#' Master function for converting prevalence rates into attack-rates
#'
#' @param X A vector of prevalence rates
#' @param PAR A set of parameters for the mechanistic model of malaria
#' @param eq Toggle for whether to calculate all attack-rates at equilibrium or
#' to only calculate the initial states at equilibrium and simulate forward
#' @param showMessages Toggle for whether the function prints messages
#' @param cpp A toggle for using the C++ version of the forward simulation
#' algorithm
#' @export
PR2AR <- function(X, PAR, eq = F, showMessages = F, cpp = T) {
    if (any(is.na(X))) {
        if (!eq) {
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
        outA <- PR2AReq(X, PAR, Bfn, showMessages)
    } else {
        outA <- simAR(X, PAR, Bfn, cpp = cpp)
        return(outA)
    }
}
