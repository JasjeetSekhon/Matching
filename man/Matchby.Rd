\name{Matchby}
\alias{Matchby}
\title{Grouped Multivariate and Propensity Score Matching}
\description{
  This function is a wrapper for the \code{\link{Match}} function which
  separates the matching problem into subgroups defined by a factor.
  This is equivalent to conducting exact matching on each level of a factor.
  Matches within each level are found as determined by the
  usual matching options.  This function is much faster for large
  datasets than the \code{\link{Match}} function itself.  For additional
  speed, consider doing matching without replacement---see the
  \code{replace} option.  This function is more limited than the
  \code{\link{Match}} function.  For example, \code{Matchby} cannot be
  used if the user wishes to provide observation specific weights.
}

\usage{
Matchby(Y, Tr, X, by, estimand = "ATT", M = 1, ties=FALSE, replace=TRUE,
        exact = NULL, caliper = NULL, AI=FALSE, Var.calc=0,
        Weight = 1, Weight.matrix = NULL, distance.tolerance = 1e-05,
        tolerance = sqrt(.Machine$double.eps), print.level=1, version="Matchby", ...)
}
\arguments{
  \item{Y}{ A vector containing the outcome of interest. Missing values are not allowed.}
  \item{Tr}{ A vector indicating the observations which are
    in the treatment regime and those which are not.  This can either be a
    logical vector or a real vector where 0 denotes control and 1 denotes
    treatment.}
  \item{X}{ A matrix containing the variables we wish to match on.
    This matrix may contain the actual observed covariates or the
    propensity score or a combination of both.}
  \item{by}{A "factor" in the sense that \code{as.factor(by)} defines the
    grouping, or a list of such factors in which case their
    interaction is used for the grouping.}
  \item{estimand}{ A character string for the estimand.  The default
    estimand is "ATT", the sample average treatment effect for the
    treated. "ATE" is the sample average treatment effect (for all), and
    "ATC" is the sample average treatment effect for the controls.}
  \item{M}{ A scalar for the number of matches which should be
    found. The default is one-to-one matching. Also see the
    \code{ties} option.}
  \item{ties}{A logical flag for whether ties should be handled
    deterministically.  By default \code{ties==TRUE}. If, for example, one
    treated observation matches more than one control observation, the
    matched dataset will include the multiple matched control observations
    and the matched data will be weighted to reflect the multiple matches.
    The sum of the weighted observations will still equal the original
    number of observations. If \code{ties==FALSE}, ties will be randomly
    broken.  \emph{If the dataset is large and there are many ties,
      setting \code{ties=FALSE} often results in a large speedup.} Whether
    two potential matches are close enough to be considered tied, is
    controlled by the \code{distance.tolerance} option.}
  \item{replace}{Whether matching should be done with replacement.  Note
    that if \code{FALSE}, the order of matches generally matters.  Matches
    will be found in the same order as the data is sorted.  Thus, the
    match(es) for the first observation will be found first and then for
    the second etc. Matching without replacement will generally increase
    bias so it is not recommended.  \emph{But if the dataset is large and
      there are many potential matches, setting \code{replace=false} often
      results in a large speedup and negligible or no bias.} Ties are
    randomly broken when \code{replace==FALSE}---see the \code{ties}
    option for details.}
  \item{exact}{ A logical scalar or vector for whether exact matching
    should be done.  If a logical scalar is provided, that logical value is
    applied to all covariates of
    \code{X}.  If a logical vector is provided, a logical value should
    be provided for each covariate in \code{X}. Using a logical vector
    allows the user to specify exact matching for some but not other
    variables.  When exact matches are not found, observations are
    dropped.  \code{distance.tolerance} determines what is considered to be an
    exact match. The \code{exact} option takes precedence over the
    \code{caliper} option.}
  \item{caliper}{ A scalar or vector denoting the caliper(s) which
    should be used when matching.  A caliper is the distance which is
    acceptable for any match.  Observations which are outside of the
    caliper are dropped. If a scalar caliper is provided, this caliper is
    used for all covariates in \code{X}.  If a vector of calipers is
    provided, a caliper value should be provide for each covariate in
    \code{X}. The caliper is interpreted to be in standardized units.  For
    example, \code{caliper=.25} means that all matches not equal to or
    within .25 standard deviations of each covariate in \code{X} are
    dropped.}
  \item{AI}{A logical flag for if the Abadie-Imbens standard error
    should be calculated. It is computationally expensive to calculate
    with large datasets. \code{Matchby} can only calculate AI SEs for ATT.
    To calculate AI errors with other estimands, please use the
    \code{\link{Match}} function.  See the \code{Var.calc} option if one
    does not want to assume homoscedasticity.}
  \item{Var.calc}{A scalar for the variance estimate
    that should be used.  By default \code{Var.calc=0} which means that
    homoscedasticity is assumed.  For values of  \code{Var.calc > 0},
    robust variances are calculated using \code{Var.calc} matches.}
  \item{Weight}{ A scalar for the type of
    weighting scheme the matching algorithm should use when weighting
    each of the covariates in \code{X}.  The default value of
    1 denotes that weights are equal to the inverse of the variances. 2
    denotes the Mahalanobis distance metric, and 3 denotes
    that the user will supply a weight matrix (\code{Weight.matrix}).  Note that
    if the user supplies a \code{Weight.matrix}, \code{Weight} will be automatically
    set to be equal to 3.}
  \item{Weight.matrix}{ This matrix denotes the weights the matching
    algorithm uses when weighting each of the covariates in \code{X}---see
    the \code{Weight} option. This square matrix should have as many
    columns as the number of columns of the \code{X} matrix. This matrix
    is usually provided by a call to the \code{\link{GenMatch}} function
    which finds the optimal weight each variable should be given so as to
    achieve balance on the covariates. \cr

    For most uses, this matrix has zeros in the off-diagonal
    cells.  This matrix can be used to weight some variables more than
    others.  For
    example, if \code{X} contains three variables and we want to
    match as best as we can on the first, the following would work well:
    \cr
    \code{> Weight.matrix <- diag(3)}\cr
    \code{> Weight.matrix[1,1] <- 1000/var(X[,1])} \cr
    \code{> Weight.matrix[2,2] <- 1/var(X[,2])} \cr
    \code{> Weight.matrix[3,3] <- 1/var(X[,3])} \cr
    This code changes the weights implied by the
    inverse of the variances by multiplying the first variable by a 1000
    so that it is highly weighted.  In order to enforce exact matching
    see the \code{exact} and \code{caliper} options.}
  \item{distance.tolerance}{This is a scalar which is used to determine if distances
    between two observations are different from zero.  Values less than
    \code{distance.tolerance} are deemed to be equal to zero.  This
    option can be used to perform a type of optimal matching}
  \item{tolerance}{ This is a scalar which is used to determine
    numerical tolerances.  This option is used by numerical routines
    such as those used to determine if a matrix is singular.}
  \item{print.level}{The level of printing. Set to '0' to turn off printing.}
  \item{version}{The version of the code to be used.  The "Matchby" C/C++
    version of the code is the fastest, and the end-user should not
    change this option.}
  \item{...}{Additional arguments passed on to \code{\link{Match}}.}
}
\details{
  \code{Matchby} is much faster for large datasets than
  \code{\link{Match}}.  But \code{Matchby} only implements a subset of
  the functionality of \code{\link{Match}}.  For example, the
  \code{restrict} option cannot be used, Abadie-Imbens standard errors
  are not provided and bias adjustment cannot be requested.
  \code{Matchby} is a wrapper for the \code{\link{Match}} function which
  separates the matching problem into subgroups defined by a factor.  This
  is the equivalent to doing exact matching on each factor, and the
  way in which matches are found within each factor is determined by the
  usual matching options. \cr

  \emph{Note that by default \code{ties=FALSE} although the default for
    the \code{Match} in \code{GenMatch} functions is \code{TRUE}.  This is
    done because randomly breaking ties in large datasets often results in
    a great speedup.}  For additional speed, consider doing matching
  without replacement which is often much faster when the dataset is
  large---see the \code{replace} option. \cr

  There will be slight differences in the matches produced by
  \code{Matchby} and \code{\link{Match}} because of how the covariates
  are weighted.  When the data is broken up into separate groups (via
  the \code{by} option), Mahalanobis distance and inverse variance
  will imply different weights than when the data is taken as whole.
}

\value{
  \item{est}{The estimated average causal effect.}
  \item{se.standard }{The usual standard error.  This is the standard error
    calculated on the matched data using the usual method of calculating
    the difference of means (between treated and control) weighted so
    that ties are taken into account.}
  \item{se }{The Abadie-Imbens standard error. This is only calculated
    if the \code{AI} option is \code{TRUE}. This standard error has
    correct coverage if \code{X} consists of either covariates or a
    known propensity score because it takes into account the uncertainty
    of the matching
    procedure.  If an estimated propensity score is used, the
    uncertainty involved in its estimation is not accounted for although the
    uncertainty of the matching procedure itself still is.}
  \item{index.treated }{A vector containing the observation numbers from
    the original dataset for the treated observations in the
    matched dataset.  This index in conjunction with \code{index.control}
    can be used to recover the matched dataset produced by
    \code{Matchby}.  For example, the \code{X} matrix used by \code{Matchby}
    can be recovered by
    \code{rbind(X[index.treated,],X[index.control,])}.}
  \item{index.control }{A vector containing the observation numbers from
    the original data for the control observations in the
    matched data.  This index in conjunction with \code{index.treated}
    can be used to recover the matched dataset produced by
    \code{Matchby}.  For example, the \code{Y} matrix for the matched dataset
    can be recovered by
    \code{c(Y[index.treated],Y[index.control])}.}
  \item{weights}{The weights for each observation in the matched
    dataset.}
  \item{orig.nobs }{The original number of observations in the dataset.}
  \item{nobs }{The number of observations in the matched dataset.}
  \item{wnobs }{The number of weighted observations in the matched dataset.}
  \item{orig.treated.nobs}{The original number of treated observations.}
  \item{ndrops}{The number of matches which were dropped because there
    were not enough observations in a given group and because of
    caliper and exact matching.}
  \item{estimand}{The estimand which was estimated.}
  \item{version}{The version of \code{\link{Match}} which was used.}
}
\references{
  Sekhon, Jasjeet S. 2011.  "Multivariate and Propensity Score
  Matching Software with Automated Balance Optimization.''
  \emph{Journal of Statistical Software} 42(7): 1-52.
  \url{https://www.jstatsoft.org/v42/i07/}

  Diamond, Alexis and Jasjeet S. Sekhon. 2013. "Genetic
  Matching for Estimating Causal Effects: A General Multivariate
  Matching Method for Achieving Balance in Observational Studies.''
  \emph{Review of Economics and Statistics}.  95 (3): 932--945.
  \url{http://sekhon.berkeley.edu/papers/GenMatch.pdf}

  Abadie, Alberto and Guido Imbens. 2006.
  ``Large Sample Properties of Matching Estimators for Average
  Treatment Effects.'' \emph{Econometrica} 74(1): 235-267.

  Imbens, Guido. 2004. Matching Software for Matlab and
  Stata.

}
\author{Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.
}
\seealso{ Also see \code{\link{Match}},
    \code{\link{summary.Matchby}},
    \code{\link{GenMatch}},
    \code{\link{MatchBalance}},
    \code{\link{balanceUV}},
    \code{\link{qqstats}}, \code{\link{ks.boot}},
    \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
}
\examples{
#
# Match exactly by racial groups and then match using the propensity score within racial groups
#

data(lalonde)

#
# Estimate the Propensity Score
#
glm1  <- glm(treat~age + I(age^2) + educ + I(educ^2) +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75, family=binomial, data=lalonde)


#save data objects
#
X  <- glm1$fitted
Y  <- lalonde$re78
Tr <- lalonde$treat

# one-to-one matching with replacement (the "M=1" option) after exactly
# matching on race using the 'by' option.  Estimating the treatment
# effect on the treated (the "estimand" option defaults to ATT).
rr  <- Matchby(Y=Y, Tr=Tr, X=X, by=lalonde$black, M=1);
summary(rr)

# Let's check the covariate balance
# 'nboots' is set to small values in the interest of speed.
# Please increase to at least 500 each for publication quality p-values.
mb  <- MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
                    hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
                    u74 + u75, data=lalonde, match.out=rr, nboots=10)

}
\keyword{nonparametric}


%  LocalWords:  MatchBalance GenMatch emph estimand ATT BiasAdjust calc dataset
%  LocalWords:  ATC ecaliper cr diag homoscedasticity rbind GerberGreenImai se
%  LocalWords:  DehejiaWahba AbadieImbens noadj cond mdata datasets wnobs url
%  LocalWords:  ndrops Abadie Imbens Econometrica Matlab Stata UC seealso Wahba
%  LocalWords:  balanceUV lalonde Dehejia psid Rajeev Sadek glm hisp
%  LocalWords:  nodegr rr nboots nmc mb Matchby ret qqstats Mahalanobis SEs
%  LocalWords:  estimands
