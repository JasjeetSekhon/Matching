test_that("Tests MatchBy ", {
suppressMessages(library(Matching))
suppressWarnings(RNGversion("3.5.3"))

data(lalonde)

X  <- cbind(lalonde$black, lalonde$age, lalonde$educ)
Y  <- lalonde$re78
Tr  <- lalonde$treat

rr2  <- Matchby(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE, by=X[,1],
                ties=TRUE, replace=TRUE, AI=TRUE)
summary(rr2)

expect_equal(rr2$est, 1907.307353)
expect_equal(rr2$se, 598.0616018)
expect_equal(all.equal(rr2$weights[1:10], c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0),
                       tolerance = 1e-6), TRUE)

rr  <- Match(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE)
summary(rr, full=TRUE)

options(digits = 20)

expect_equal(rr$est[1,1], 1907.307353)
expect_equal(rr$se, 598.0616018)
expect_equal(all.equal(rr$weights[1:10], c(1.0000000, 0.3333333, 0.3333333,
                                           0.3333333, 0.3333333, 0.3333333,
                                           0.3333333, 1.0000000, 0.3333333, 0.3333333),
                       tolerance = 1e-6), TRUE)

expect_equal((abs(rr$est-rr2$est) < 1e-10)[1,1], TRUE)
expect_equal((abs(rr$se-rr2$se) < 1e-10), TRUE)
expect_equal((abs(rr$se.standard-rr2$se.standard) < 1e-10), TRUE)

})
