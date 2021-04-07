\name{qqstats}
\alias{qqstats}
\title{QQ Summary Statistics}
\description{
  This function calculates a set of summary statistics for the QQ
  plot of two samples of data.  The summaries are useful for determining
  if the two samples are from the same distribution.  If
  \code{standardize==TRUE}, the empirical CDF is used instead of the
  empirical-QQ plot.  The later retains the scale of the variable.
}
\usage{
qqstats(x, y, standardize=TRUE, summary.func)
}
\arguments{
  \item{x}{The first sample.}
  \item{y}{The second sample.}
  \item{standardize}{A logical flag for whether the statistics should be
  standardized by the empirical cumulative distribution functions of the
  two samples.}
  \item{summary.func}{A user provided function to summarize the
    difference between the two distributions.  The function should
    expect a vector of the differences as an argument and return summary
    statistic.  For example, the \code{\link{quantile}} function is a
    legal function to pass in.}
}
\value{
  \item{meandiff}{The mean difference between the QQ plots of the two
    samples.}
  \item{mediandiff}{The median difference between the QQ plots of the
    two samples.}
  \item{maxdiff}{The maximum difference between the QQ plots of the two
    samples.}
  \item{summarydiff}{If the user provides a \code{summary.func}, the
    user requested summary difference is returned.}
  \item{summary.func}{If the user provides a \code{summary.func}, the
    function is returned.}
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.
}
\references{
  Sekhon, Jasjeet S. 2011.  "Multivariate and Propensity Score
  Matching Software with Automated Balance Optimization.''
  \emph{Journal of Statistical Software} 42(7): 1-52.
  \url{https://www.jstatsoft.org/v42/i07/}

  Diamond, Alexis and Jasjeet S. Sekhon. Forthcoming. "Genetic Matching for
  Estimating Causal Effects: A General Multivariate Matching Method for
  Achieving Balance in Observational Studies.'' \emph{Review of Economics and Statistics}.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}
  }
  \seealso{ Also see \code{\link{ks.boot}},
    \code{\link{balanceUV}}, \code{\link{Match}},
    \code{\link{GenMatch}},
    \code{\link{MatchBalance}},
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
qqout  <- qqstats(lalonde$re75[rr$index.treated], lalonde$re75[rr$index.control])
print(qqout)
}
\keyword{htest}
\keyword{distribution}

%  LocalWords:  Kolmogorov nboots ksboot pvalue UC url Abadie emph seealso psid
%  LocalWords:  balanceUV GenMatch MatchBalance GerberGreenImai Wahba
%  LocalWords:  lalonde Dehejia Rajeev Sadek glm hisp nodegr estimand rr htest
%  LocalWords:  maxdiff QQ summarydiff func qqout qqstats meandiff mediandiff
