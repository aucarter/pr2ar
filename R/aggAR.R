#' Aggregates attack-rates from time step to annual
#'
#' @param A A vector of attack-rates
#' @param dt length of time step
#' @export
aggAR <- function(A, dt) {
    len = round(365 / dt)
    Aout = c()
    for(i in 1:floor(length(A) / len)) {
        aggA = 1 - prod(1 - A[((i - 1) * len + 1): (i * len)])
        Aout = c(Aout, aggA)
    }
    return(Aout)
}
