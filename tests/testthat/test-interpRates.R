test_that("multiplication works", {
    inputRates <- seq(0.5, 0.9, 0.1)
    inputStep = 365
    outputStep = 14
    expectedN  = (length(inputRates) - 1) * round(inputStep / outputStep) + 1
    interpN = length(interpRates(rates = inputRates,
                                 inStep = inputStep,
                                 outStep = outputStep))
    expect_equal(expectedN, interpN)
})
