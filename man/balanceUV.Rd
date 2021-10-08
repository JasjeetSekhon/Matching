\name{balanceUV}
\alias{balanceUV}
\title{Univariate Balance Tests}
\description{
  This function provides a number of univariate balance metrics.
  Generally, users should call \code{\link{MatchBalance}} and not this function
  directly.
}
\usage{
balanceUV(Tr, Co, weights = rep(1, length(Co)), exact = FALSE, ks=FALSE,
          nboots = 1000, paired=TRUE, match=FALSE,
          weights.Tr=rep(1,length(Tr)), weights.Co=rep(1,length(Co)),
          estimand="ATT")
}
\arguments{
  \item{Tr}{A vector containing the treatment observations.}
  \item{Co}{A vector containing the control observations.}
  \item{weights}{A vector containing the observation specific
    weights. Only use this option when the treatment and control
    observations are paired (as they are after matching).}
  \item{exact}{A logical flag indicating if the exact Wilcoxon test
    should be used instead of the test with a correction.  See
    \code{\link{wilcox.test}} for details.}
  \item{ks}{ A logical flag for if the univariate bootstrap
    Kolmogorov-Smirnov (KS) test should be calculated.  If the ks option
    is set to true, the univariate KS test is calculated for all
    non-dichotomous variables.  The bootstrap KS test is consistent even
    for non-continuous variables.  See \code{\link{ks.boot}} for more
    details.}
  \item{nboots}{The number of bootstrap samples to be run for the
    \code{ks} test.  If zero, no bootstraps are done.  Bootstrapping is
    highly recommended because the bootstrapped Kolmogorov-Smirnov test
    only provides correct coverage even for non-continuous covariates. At
    least 500 \code{nboots} (preferably 1000) are
    recommended for publication quality p-values.}
  \item{paired}{A flag for if the paired \code{\link{t.test}} should be used.}
  \item{match}{A flag for if the \code{Tr} and \code{Co} objects are the result
    of a call to \code{\link{Match}}.}
  \item{weights.Tr}{A vector of weights for the treated observations.}
  \item{weights.Co}{A vector of weights for the control observations.}
  \item{estimand}{This determines if the standardized mean difference
    returned by the \code{sdiff} object is standardized by the variance of
    the treatment observations (which is done if the estimand is either
    "ATE" or "ATT") or by the variance of the control observations (which
    is done if the estimand is "ATC").}
}
\value{
  \item{sdiff}{This is the standardized difference between the treated
    and control units multiplied by 100.  That is, 100 times the mean
    difference  between treatment and control units divided by the standard
    deviation of the treatment
    observations alone if the estimand is either \code{ATT} or
    \code{ATE}.  The variance of the control observations are used if
    the estimand is \code{ATC}.}
  \item{sdiff.pooled}{This is the standardized difference between the
    treated and control units multiplied by 100 using the pooled
    variance.  That is, 100 times the mean
    difference between treatment and control units divided by the pooled standard deviation as in
    Rosenbaum and Rubin (1985).}
  \item{mean.Tr}{The mean of the treatment group.}
  \item{mean.Co}{The mean of the control group.}
  \item{var.Tr}{The variance of the treatment group.}
  \item{var.Co}{The variance of the control group.}
  \item{p.value}{The p-value from the two-sided weighted \code{\link{t.test}}.}
  \item{var.ratio}{var.Tr/var.Co.}
  \item{ks}{The object returned by \code{\link{ks.boot}}.}
  \item{tt}{The object returned by two-sided weighted
    \code{\link{t.test}}.}
  \item{qqsummary}{The return object from a call to
    \code{\link{qqstats}} with standardization---i.e., balance test
    based on the empirical CDF.}
  \item{qqsummary.raw}{The return object from a call to
    \code{\link{qqstats}} without standardization--i.e., balance tests
    based on the empirical QQ-plot which retain the scale of the
    variable.}
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.
}
\references{
  Sekhon, Jasjeet S. 2011.  "Multivariate and Propensity Score
  Matching Software with Automated Balance Optimization.''
  \emph{Journal of Statistical Software} 42(7): 1-52.
  \doi{10.18637/jss.v042.i07}

  Diamond, Alexis and Jasjeet S. Sekhon. 2013. "Genetic
  Matching for Estimating Causal Effects: A General Multivariate
  Matching Method for Achieving Balance in Observational Studies.''
  \emph{Review of Economics and Statistics}.  95 (3): 932--945.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}

  Rosenbaum, Paul R. and Donald B. Rubin. 1985. ``Constructing a Control
  Group Using Multivariate Matched Sampling Methods That Incorporate the
  Propensity Score.'' \emph{The American Statistician} 39:1 33-38.

  Hollander, Myles and Douglas A. Wolfe. 1973. \emph{Nonparametric
    statistical inference}. New York: John Wiley & Sons.
  }
  \seealso{ Also see \code{\link{summary.balanceUV}}, \code{\link{qqstats}}
    \code{\link{ks.boot}}, \code{\link{Match}}, \code{\link{GenMatch}},
    \code{\link{MatchBalance}},
    \code{\link{GerberGreenImai}}, \code{\link{lalonde}} }
\examples{
data(lalonde)
attach(lalonde)

foo  <- balanceUV(re75[treat==1],re75[treat!=1])
summary(foo)

}
\keyword{univar}

%  LocalWords:  balanceUV Univariate MatchBalance nboots Wilcoxon wilcox sdiff
%  LocalWords:  Kolmogorov tt UC url emph seealso GenMatch lalonde
%  LocalWords:  GerberGreenImai univar qqstats
