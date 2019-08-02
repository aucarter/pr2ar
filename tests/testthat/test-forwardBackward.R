test_that("Prevalence generated from a specific attack-rate can be recovered", {
    n = 20
    Bfn = makeBdrugs_age
    PAR = list(In = 5, Cn = 2, Q = 0.95)
    # Make up attack-rate and treatment rate
    AR <- (1 + 0.1*sin(seq(0, 10, length.out = n + 1)))*seq(0.05, 0.4, length.out = n + 1)
    Tx <- seq(0.05, 0.3, length.out = n)
    # Generate a corresponding vector of prevalence
    D <- makeD(PAR)
    PAR$A <- AR[1]
    PAR$rho <- Tx[1]
    Y <- findYeq(PAR, Bfn)
    X <- c()
    for(i in 1:length(Tx)) {
        PAR$A = AR[i + 1]
        PAR$rho = Tx[i]
        B = Bfn(PAR)
        Y = B%*%Y
        X = c(X, D%*%Y)
    }
    AR2 = PR2AR(X, Tx, PAR, cpp = T)
    AR3 = PR2AR(X, Tx, PAR, cpp = F)
    expect_lt(sum(head(AR[-1], -1) - AR2), 0.0001, label = "C++ version")
    expect_lt(sum(head(AR[-1], -1) - AR3), 0.0001, label = "R version")
})
