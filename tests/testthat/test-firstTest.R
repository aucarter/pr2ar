test_that("models generate Markov matrices", {
    PAR <- list(A=.3, Q=0.95, rho=0.05, In=5, Cn=2, d = 0)
    B = makeBdrugs_age(PAR)
    expect_equal(sum(colSums(B) - 1) < .Machine$double.eps, T)
})
