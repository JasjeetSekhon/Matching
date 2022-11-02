test_that("Tests GenMatch", {
suppressMessages(library(rgenoud))
suppressMessages(library(Matching))
suppressWarnings(RNGversion("3.5.3"))

set.seed(3101)
data(lalonde)
attach(lalonde)

#The covariates we want to match on
X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74)

#The covariates we want to obtain balance on
BalanceMat <- cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74,
                    I(re74*re75))

#
#Let's call GenMatch() to find the optimal weight to give each
#covariate in 'X' so as we have achieved balance on the covariates in
#'BalanceMat'. This is only an example so we want GenMatch to be quick
#so the population size has been set to be only 16 via the 'pop.size'
#option. This is *WAY* too small for actual problems.
#For details see http://sekhon.berkeley.edu/papers/MatchingJSS.pdf.
#
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATE", M=1,
                   pop.size=16, max.generations=10, wait.generations=1,
                   unif.seed=3392, int.seed=8282, print.level=0)
print(genout)

#The outcome variable
Y=re78/1000

#
# Now that GenMatch() has found the optimal weights, let's estimate
# our causal effect of interest using those weights
#
mout <- Match(Y=Y, Tr=treat, X=X, estimand="ATE", Weight.matrix=genout)
summary(mout)

#
#Let's determine if balance has actually been obtained on the variables of interest
#
mb <- MatchBalance(treat~age +educ+black+ hisp+ married+ nodegr+ u74+ u75+
                     re75+ re74+ I(re74*re75),
                   match.out=mout, nboots=500)

expect_equal(nrow(genout$Weight.matrix), 10)
expect_equal(ncol(genout$Weight.matrix), 10)

# For more examples see: http://sekhon.berkeley.edu/matching/R.
before_sdiff <- unlist(lapply(mb$BeforeMatching, function(x){return(x$sdiff)}))
after_sdfiff <- unlist(lapply(mb$AfterMatching, function(x){return(x$sdiff)}))

results <- data.frame(before = before_sdiff,
                      after = after_sdfiff)

results_saved <- data.frame(before = c(10.655038555334192,  12.806026640590549,   4.476700682905446,
                                       -20.340738543987911,   8.999512178599344, -27.750943560007322,
                                       -9.189507193658958, -17.225298606861006, 8.236275506237831,
                                       -0.234370781139797,  -2.779905583506529),
                            after = c(0.698105986659637, -1.579765862957684, -0.609524006759037,
                                      0.000000000000000,  0.000000000000000,  0.000000000000000,
                                      -0.507141556044142, -0.939497016216188,  2.263685815459904,
                                      -0.280105002856819, -2.639651943919278))

# Skip the exact tests if not on MacOS
skip_if_not_mac()

# Check the sdiff before matching
expect_equal(all.equal(results[,1], results_saved[,1]), TRUE)
# Check the sdiff after matching
expect_equal(all.equal(results[,2], results_saved[,2]), TRUE)


})
