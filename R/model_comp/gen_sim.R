#' Simulate a time-series of attack-rates
#'
#' @param type A string with the type of simulation to generate
#' @param step Time step in days
#' @param years Length of simulation in years
#' @param plot Boolean for including a plot of the simulated AR
#' @example genPR(PAR= list(In = 5, Cn = 2, Q = 0.95, dt = 10, rho = 0.5),
#'          Bfn = makeBdrugs_age,
#'          A_type = "seasonal",
#'          Tx_type = "increasing",
#'          plot = T)
#'
#' @export
genA <- function(type, step = 10, years = 10, plot = F) {
    if(type == "seasonal") {
        level <- 0.3
        amp <- 0.1
        wavelen <- 0.5 # per year
        time <- seq(0, 365 * years - (365 * years) %% step,  step)
        A <- level + amp * sin(2 * pi * time / (365 * wavelen))
    }
    if(plot) {
        plot(time / 360, A, type = "l")
    }
    return(A)
}
#' Simulate a time-series of treatment rates
#'
#' @param type A string with the type of simulation to generate
#' @param step Time step in days
#' @param years Length of simulation in years
#' @param plot Boolean for including a plot of the simulated Tx
#'
#' @export
genTx <- function(type, step = 10, years = 10, plot = F) {
    if(type == "increasing") {
        min <- 0.1
        max <- 0.4
        time <- c(0, 365 * years - (365 * years) %% step)
        Tx <- approx(x = time, y = c(min, max), n = ((365 * years - (365 * years) %% step) / step + 1))$y
    }
    if(plot) {
        plot(time / 360, Tx, type = "l")
    }
    return(Tx)
}

#' Simulate aprevalence rates with a given model and set of attack-rates and treatment rates
#'
#' @param PAR A set of model parameters
#' @param Bfn A function for building the matrix represention of the system of
#' equations that update the state vector
#' @param A_type A string with the type of attack-rate simulation to generate
#' @param Tx_type A string with the type of treatment rate simulation to generate
#' @param step Time step in days
#' @param years Length of simulation in years
#' @param plot Boolean for including a plot of the simulated PR
#'
#' @export

genPR<- function(PAR, Bfn, A_type, Tx_type, step = 10, years = 10, plot = F) {
    A <- genA(A_type, step, years)
    Tx <- genTx(Tx_type, step, years)
    PR <- AR2PR(A, Tx, PAR, Bfn)$X
    if(plot) {
        time <- seq(0, 365 * years - (365 * years) %% step,  step)
        plot(time / 360, PR, type = "l")
    }
    return(PR)
}
