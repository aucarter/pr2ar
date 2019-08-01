test_that("Prevalence generated from a specific attack-rate can be recovered", {
    Bfn = makeBdrugs_age
    PAR = list(rho = 0.1, In = 5, Cn = 2, Q = 0.95)
    # Make up attack-rate
    AR <- (1 + 0.1*sin(seq(0, 10, length.out = 20)))*seq(0.05, 0.4, length.out = 20)
    # Generate a corresponding vector of prevalence
    D <- makeD(PAR)
    PAR$A <- AR[1]
    Y <- findYeq(PAR, Bfn)
    X <- c(D%*%Y)
    for(i in 2:length(AR)) {
        PAR$A = AR[i]
        B = Bfn(PAR)
        Y = B%*%Y
        X = c(X, D%*%Y)
    }
    AR2 = PR2AR(X, PAR)$A
    expect_less_than(sum(head(AR, -1) - AR2), 0.0001)
})
