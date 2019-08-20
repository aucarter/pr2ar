#' Master function for converting prevalence rates into attack-rates
#'
#' @param X A vector of prevalence rates
#' @param Tx A vector of case-management rates
#' @param PAR A set of parameters for the mechanistic model of malaria
#' @param Xinterval An integer with the number of days between PR observations
#' @param eq Toggle for whether to calculate all attack-rates at equilibrium or
#' to only calculate the initial states at equilibrium and simulate forward
#' @param showMessages Toggle for whether the function prints messages
#' @param cpp A toggle for using the C++ version of the forward simulation
#' algorithm
#' @export
PR2AR <- function(X, Tx, PAR, Xinterval = 365, eq = F, showMessages = F, cpp = T) {
    # Screen for NAs
    if (any(is.na(X))) {
        if (!eq) {
            return(rep(NA, length(X)))
        } else {
            return(list(A = rep(NA, length(X)), Y = matrix(NA, nrow = length(makeD(PAR)), ncol = length(X))))
        }
    }

    # Interpolate inputs to match time step
    X = interpRates(X, inStep = Xinterval, outStep = PAR$dt)
    Tx = interpRates(Tx, inStep = Xinterval, outStep = PAR$dt)

    # Choose model
    if (PAR$In > 1 & PAR$Cn > 0) {
        Bfn = makeBdrugs_age
    } else if (PAR$In > 1) {
        Bfn = makeBage
    } else if (PAR$Cn > 0) {
        Bfn = makeBdrugs
    } else {
        Bfn = makeB
    }

    # Calculate AR
    if (eq) {
        AR <- PR2AReq(X, Tx, PAR, Bfn, showMessages)
    } else {
        AR <- simAR(X, Tx, PAR, Bfn, cpp = cpp)
    }

    # Return to original input time step
    outA = aggAR(AR, inStep = Xinterval, outStep = PAR$dt)
    return(outA)
}
