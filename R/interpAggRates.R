#' Interpolate a vector of rates to match the model time step
#'
#' @param rates A vector of rates
#' @param inStep Input interval size in days
#' @param outStep output interval size in days
#' @export
interpRates <- function(rates, inStep = 365, outStep = 10) {
    size = (length(rates) - 1) * round(inStep / outStep) + 1
    if(!any(is.na(rates))) {
        smooth <- stats::smooth.spline(x = 1:length(rates), y = rates)
        interpX = seq(1, length(rates), length.out =  size)
        interpR <- stats::predict(smooth$fit, interpX)$y
        interpR[interpR < 0] = 0
    } else {
        interpR <- rep(NA, size)
    }
    return(interpR)
}

#' Aggregates attack-rates from time step to annual
#'
#' @param A A vector of attack-rates
#' @param inStep Input interval size in days
#' @param outStep output interval size in days
#' @export
aggAR <- function(A, inStep, outStep) {
    len = round(inStep / outStep)
    Aout = c()
    for(i in 1:floor(length(A) / len)) {
        aggA = 1 - prod(1 - A[((i - 1) * len + 1):(i * len)])
        Aout = c(Aout, aggA)
    }
    return(Aout)
}
