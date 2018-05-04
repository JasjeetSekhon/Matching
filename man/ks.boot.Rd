\name{ks.boot}
\alias{ks.boot}
\title{Bootstrap Kolmogorov-Smirnov}
\description{
  This function executes a bootstrap version of the univariate
  Kolmogorov-Smirnov test which provides correct coverage even when
  the distributions being compared are not entirely continuous.  Ties
  are allowed with this test unlike the traditional Kolmogorov-Smirnov
  test.
}
\usage{
ks.boot(Tr, Co, nboots=1000, alternative = c("two.sided","less","greater"),
        print.level=0)
}
\arguments{
  \item{Tr}{A vector containing the treatment observations.}
  \item{Co}{A vector containing the control observations.}
  \item{nboots}{The number of bootstraps to be performed.  These are, in
    fact, really Monte Carlo simulations which are preformed in order to
    determine the proper p-value from the empiric.}
  \item{alternative}{indicates the alternative hypothesis and must be one of
    '"two.sided"' (default), '"less"', or '"greater"'.  You can
    specify just the initial letter.  See \code{\link{ks.test}} for details.}
  \item{print.level}{If this is greater than 1, then the simulation count is
    printed out while the simulations are being done.}
}
\value{
  \item{ks.boot.pvalue}{The bootstrap p-value of the Kolmogorov-Smirnov test
    for the hypothesis that the probability densities for both the
    treated and control groups are the same.}
  \item{ks}{Return object from \code{\link{ks.test}}.}
  \item{nboots}{The number of bootstraps which were completed.}
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.  
}
\references{
  Sekhon, Jasjeet S. 2011.  "Multivariate and Propensity Score
  Matching Software with Automated Balance Optimization.''
  \emph{Journal of Statistical Software} 42(7): 1-52.
  \url{http://www.jstatsoft.org/v42/i07/}

  Diamond, Alexis and Jasjeet S. Sekhon. 2013. "Genetic Matching for
  Estimating Causal Effects: A General Multivariate Matching Method for
  Achieving Balance in Observational Studies.'' \emph{Review of
  Economics and Statistics}.  95 (3): 932--945.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}
  
  Abadie, Alberto.  2002.  ``Bootstrap Tests for Distributional Treatment
  Effects in Instrumental Variable Models.'' \emph{Journal of the
  American Statistical Association}, 97:457 (March) 284-292.
  }
  \seealso{ Also see \code{\link{summary.ks.boot}},
    \code{\link{qqstats}}, \code{\link{balanceUV}}, \code{\link{Match}},
    \code{\link{GenMatch}}, \code{\link{MatchBalance}}, 
    \code{\link{GerberGreenImai}}, \code{\link{lalonde}} }
\examples{
#
# Replication of Dehejia and Wahba psid3 model
#
# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in
# Non-Experimental Studies: Re-Evaluating the Evaluation of Training
# Programs.''Journal of the American Statistical Association 94 (448):
# 1053-1062.
#
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
# Estimating the treatment effect on the treated (the "estimand" option which defaults to 0).
#
rr  <- Match(Y=Y,Tr=Tr,X=X,M=1);
summary(rr)

#
# Do we have balance on 1975 income after matching?
#
ks  <- ks.boot(lalonde$re75[rr$index.treated], lalonde$re75[rr$index.control], nboots=500)
summary(ks)
}
\keyword{htest}

%  LocalWords:  Kolmogorov nboots ksboot pvalue UC url Abadie emph seealso psid
%  LocalWords:  balanceUV GenMatch MatchBalance GerberGreenImai Wahba
%  LocalWords:  lalonde Dehejia Rajeev Sadek glm hisp nodegr estimand rr htest
