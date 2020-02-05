suppressMessages(library(Matching))
suppressWarnings(RNGversion("3.5.3"))

data(lalonde)

X  <- cbind(lalonde$black, lalonde$age, lalonde$educ)
Y  <- lalonde$re78
Tr  <- lalonde$treat

rr2  <- Matchby(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE, by=X[,1],
                ties=TRUE, replace=TRUE, AI=TRUE)
summary(rr2)

rr  <- Match(Y=Y, Tr=Tr, X=X, M=1, exact=TRUE)
summary(rr, full=TRUE)

abs(rr$est-rr2$est) < 1e-10
abs(rr$se-rr2$se) < 1e-10
abs(rr$se.standard-rr2$se.standard) < 1e-10
