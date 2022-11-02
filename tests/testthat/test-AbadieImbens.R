test_that("Tests replication of Guido Imbens lalonde_exper_04feb2.m file", {
# Replication of Guido Imbens lalonde_exper_04feb2.m file
# See http://elsa.berkeley.edu/~imbens/estimators.shtml
# with balance checks

suppressMessages(library(Matching))
suppressWarnings(RNGversion("3.5.3"))

data(lalonde)

X  <- lalonde$age
Z  <- X;
V  <- lalonde$educ;
Y  <- lalonde$re78/1000;
T  <- lalonde$treat;
w.educ=exp((lalonde$educ-10.1)/2);

res  <- matrix(nrow=1,ncol=3)

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)

expect_equal(rr$est[1,1], 1.7851911148121248907)
expect_equal(rr$se, 0.6867166773225281684)
expect_equal(all.equal(rr$weights[1:10], c(0.14285714285714284921, 0.14285714285714284921, 0.14285714285714284921,
                                           0.14285714285714284921, 0.14285714285714284921, 0.14285714285714284921,
                                           0.14285714285714284921, 0.06250000000000000000, 0.06250000000000000000, 0.06250000000000000000)),TRUE)

res[1,]  <- cbind(1,rr$est,rr$se)


X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74, lalonde$re75)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)

res  <- rbind(res,cbind(2,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=3,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(4,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATT",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(5,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATC",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(6,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=2,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(7,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=3,Var.calc=0,
             Weight.matrix=diag(4), sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(8,rr$est,rr$se))


rr  <- Match(Y=Y,Tr=T,X=X,Z=X,V=V,estimand="ATE",M=1,BiasAdj=TRUE,Weight=1,Var.calc=0,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(9,rr$est,rr$se))

Z  <- cbind(lalonde$married, lalonde$age)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=TRUE,Weight=1,Var.calc=0,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(10,rr$est,rr$se))

V  <- lalonde$age
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,exact=TRUE,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(11,rr$est,rr$se))

V  <- cbind(lalonde$married, lalonde$u74)
rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,exact=TRUE,
             sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(12,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,sample=FALSE);
summary(rr)
res  <- rbind(res,cbind(13,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=3,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(14,rr$est,rr$se))

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,estimand="ATE",M=1,BiasAdj=FALSE,Weight=1,Var.calc=0,
             weights=w.educ,sample=TRUE);
summary(rr)
res  <- rbind(res,cbind(15,rr$est,rr$se))


V  <- lalonde$age
Z  <- cbind(lalonde$married, lalonde$age)
X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74, lalonde$re75)
weight  <- w.educ
Weight.matrix  <- diag(4)

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,
             sample=FALSE, M=3, estimand="ATT", BiasAdj=TRUE, Weight=3, exact=TRUE,Var.calc=3,
             weights=w.educ, Weight.matrix=Weight.matrix);
summary(rr)
res  <- rbind(res,cbind(75,rr$est,rr$se))


V  <- lalonde$married;
Z  <- cbind(lalonde$age, lalonde$re75);
X  <- cbind(lalonde$age, lalonde$educ, lalonde$re74);

rr  <- Match(Y=Y,Tr=T,X=X,Z=Z,V=V,
             sample=TRUE, M=3, estimand="ATE", BiasAdj=TRUE, Weight=2, exact=TRUE,Var.calc=0,
             weights=w.educ);
summary(rr)
res  <- rbind(res,cbind(76,rr$est,rr$se))

cat("\nResults:\n")
#print(res)

res_test <- data.frame(Num = c(1:16),
                       est = c(1.78519111481212, 1.71440714918138, 1.53630268121115, 1.72691301785071,
                               1.70550874262821, 1.59261673326378, 1.71440714918138, 1.63091519500783,
                               1.71968656670868, 1.36928912506340, 1.36928912506340, 1.71440714918138,
                               1.71440714918138, 3.14695956932726, 3.16566014738572, 3.50457962724634),
                       se = c(0.686716677322528, 0.740129179632503, 0.661930423623259, 0.836629632308480,
                              0.820095935752824, 0.684729829415588, 0.740129179632503, 0.745233606420550,
                              0.743608513017050, 0.404228086032945, 0.404228086032945, 0.743384632957198,
                              0.696272433717223, 1.070914907138643, 0.265657268489384, 0.574613885901435))

# Check results with testthat
for ( i in 1:nrow(res)) {
  # Check est against saved versions
  expect_equal(res[i,2], res_test[i,2])
  # Check SE against saved versions
  expect_equal(res[i,3], res_test[i,3])
}

})
