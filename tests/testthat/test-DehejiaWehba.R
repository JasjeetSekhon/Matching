test_that("Tests replication of Dehejia and Wahba psid3 model", {
#
# Replication of Dehejia and Wahba psid3 model
#
# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in Non-Experimental Studies: Re-Evaluating the
# Evaluation of Training Programs.''Journal of the American Statistical Association 94 (448): 1053-1062.
#

suppressMessages(library(Matching))
suppressWarnings(RNGversion("3.5.3"))

# Replication of Dehejia and Wahba psid3 model.

# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in
# Non-Experimental Studies: Re-Evaluating the # Evaluation of Training
# Programs.''Journal of the American Statistical Association 94 (448):
# 1053-1062.

set.seed(10391)
data(lalonde)

#
# Estimate the propensity model
#
glm1  <- glm(treat~age + I(age^2) + educ + I(educ^2) + black +
               hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
               u74 + u75, family=binomial, data=lalonde)

#
#save data objects
#
X  <- glm1$fitted
Y  <- lalonde$re78
Tr  <- lalonde$treat

#
# one-to-one matching with replacement (the "M=1" option).
# Estimating the treatment effect on the treated (the "estimand" option which defaults ATT).
#
rr  <- Match(Y=Y,Tr=Tr,X=X,M=1);
summary(rr)

#
# Let's check for balance
#
mb  <- MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
                      hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
                      u74 + u75, data=lalonde, match.out=rr, nboots=0)

before_sdiff <- unlist(lapply(mb$BeforeMatching, function(x){return(x$sdiff)}))
after_sdfiff <- unlist(lapply(mb$AfterMatching, function(x){return(x$sdiff)}))

results <- data.frame(before = before_sdiff,
                      after = after_sdfiff)

results_saved <- data.frame(before = c(10.655038555334192,   9.293693091477756,  12.806026640590549,
                                       17.012014396718083,   4.476700682905446, -20.340738543987911,
                                       8.999512178599344, -27.750943560007322, -0.234370781139797,
                                       -7.472147883152990,   8.236275506237831,   2.602407688369587,
                                       -9.189507193658958, -17.225298606861006),
                      after = c(11.31676699418897, 10.27537464950881, -6.67488011805762,
                                -5.46597309818458, -4.44818666581684,  4.55913105296281,
                                5.73498325106821,  3.55722859109378,  9.64392627946963,
                                13.16663153873934,  7.28274018652368,  6.70756470721747,
                                5.16080465438050, -4.21815256629930))

for ( i in 1:nrow(results)) {
  # Check the sdiff before matching
  expect_equal(results[i,1], results_saved[i,1])
  # Check the sdiff after matching
  expect_equal(results[i,2], results_saved[i,2])
}


})
