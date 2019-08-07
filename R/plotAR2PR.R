#' Make plots showing the relationship between AR and PR for a number of different Tx levels
#'
#' @param PAR Model parameters
#' @param Bfn Function for the model
#' @param AR_range Min and max of the range of AR plotted
#' @param logX Toggle for taking the log of the x-axis
#' @param logY Toggle for taking the log of the y-axis
#' @param plotC Toggle for plotting the proportion treated
#' @param plotXover1_C Toggle for plotting PR relative to those not on Tx
#' @export
plotAR2PR <- function(PAR, Bfn, AR_range = c(0.01, 1), logX = T, logY = F,
                      plotC = T, plotXover1_C = F) {
    # Prep data
    AR = seq(AR_range[1], AR_range[2], length.out = 100)
    rho = c(seq(0, 0.075, 0.025), seq(0.1, 0.8, .1))
    dt <- data.table()
    for(rho_i in rho) {
        PAR$rho <- rho_i
        PR <- AR2PR(AR, Tx = rep(rho_i, length(AR)), PAR, Bfn)
        X_i = PR$X
        C_i = PR$C
        add.dt <-  data.table(AR = AR, rho = rho_i, X = X_i, C = C_i)
        dt <- rbind(dt, add.dt)
    }

    # Plot
    AR <- X <- C <- NULL # annoying thing for checker
     gg <- ggplot(dt, aes(x = AR, y = X)) + geom_line() + facet_wrap(~rho) +
         theme_classic() +
         labs(title = "Attack-rate and PfPR for different case management rates",
              subtitle = paste0("Anti-malarial drug use w/o confirmed malaria: ", PAR$d)) +
         xlab("Attack-rate") + ylab("PfPR")

    if(logX) {
        gg = gg + scale_x_continuous(trans='log10')
    }
    if(logY) {
        gg = gg + scale_y_continuous(trans='log10')
    } else {
        gg = gg + ylim(c(0, 1))
    }
    if(plotC) {
        top.dt <- copy(dt)[, C := 0]
        top.dt <- top.dt[nrow(top.dt) - 1:nrow(top.dt) + 1]
        poly.dt <- rbind(dt, top.dt)
        gg <- gg + geom_polygon(data = poly.dt, aes(x= AR, y = (1-C), fill = "Treated"), alpha = 0.5) +
             theme(legend.position="bottom", legend.title = element_blank())
    }
    if(plotXover1_C) {
        gg <- gg + geom_line(aes(x = AR, y = X / (1 - C)), color = "red")
    }
    gg
}
