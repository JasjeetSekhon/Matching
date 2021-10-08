\name{Match}
\alias{Match}
\title{Multivariate and Propensity Score Matching Estimator for Causal Inference}
\description{
  \code{Match} implements a variety of algorithms for multivariate
  matching including propensity score, Mahalanobis and inverse variance
  matching.  The function is intended to be used in conjunction with the
  \code{MatchBalance} function which determines the extent to which
  \code{Match} has been able to achieve covariate balance.  In order to
  do propensity score matching, one should estimate the propensity model
  before calling \code{Match}, and then send \code{Match} the propensity
  score to use.  \code{Match} enables a wide variety of matching
  options including matching with or without replacement, bias
  adjustment, different methods for handling ties, exact and caliper
  matching, and a method for the user to fine tune the matches via a
  general restriction matrix.  Variance estimators include the usual
  Neyman standard errors, Abadie-Imbens standard errors, and robust
  variances which do not assume a homogeneous causal effect. The
  \code{\link{GenMatch}} function can be used to \emph{automatically
  find balance} via a genetic search algorithm which determines the
  optimal weight to give each covariate.
}
\usage{
Match(Y=NULL, Tr, X, Z = X, V = rep(1, length(Y)), estimand = "ATT", M = 1,
      BiasAdjust = FALSE, exact = NULL, caliper = NULL, replace=TRUE, ties=TRUE,
      CommonSupport=FALSE,Weight = 1, Weight.matrix = NULL, weights = NULL,
      Var.calc = 0, sample = FALSE, restrict=NULL, match.out = NULL,
      distance.tolerance = 1e-05, tolerance=sqrt(.Machine$double.eps),
      version="standard")
}
\arguments{
  \item{Y}{ A vector containing the outcome of interest.
    Missing values are not allowed.  An outcome vector is not required
    because the matches generated will be the same regardless of the
    outcomes.  Of course, without any outcomes no causal effect
    estimates will be produced, only a matched dataset. }
  \item{Tr}{ A vector indicating the observations which are
    in the treatment regime and those which are not.  This can either be a
    logical vector or a real vector where 0 denotes control and 1 denotes
    treatment.}
  \item{X}{ A matrix containing the variables we wish to match on.
    This matrix may contain the actual observed covariates or the
    propensity score or a combination of both. All columns of this
    matrix must have positive variance or \code{Match} will return an
    error.}
  \item{Z}{ A matrix containing the covariates for which we wish to make
    bias adjustments.}
  \item{V}{ A matrix containing the covariates for which the variance
    of the causal effect may vary. Also see the \code{Var.calc} option,
    which takes precedence.}
  \item{estimand}{ A character string for the estimand.  The default
    estimand is "ATT", the sample average treatment effect for the
    treated. "ATE" is the sample average treatment effect, and
    "ATC" is the sample average treatment effect for the controls.}
  \item{M}{A scalar for the number of matches which should be
    found. The default is one-to-one matching. Also see the \code{ties}
    option.}
  \item{BiasAdjust}{ A logical scalar for whether regression adjustment
    should be used. See the \code{Z} matrix.}
  \item{exact}{ A logical scalar or vector for whether exact matching
    should be done. If a logical scalar is provided, that logical value is
      applied to all covariates in
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
    provided, a caliper value should be provided for each covariate in
    \code{X}. The caliper is interpreted to be in standardized units.  For
    example, \code{caliper=.25} means that all matches not equal to or
    within .25 standard deviations of each covariate in \code{X} are
    dropped. Note that dropping observations generally changes the
    quantity being estimated.}
  \item{replace}{A logical flag for whether matching should be done with
    replacement.  Note that if \code{FALSE}, the order of matches
    generally matters.  Matches will be found in the same order as the
    data are sorted.  Thus, the match(es) for the first observation will
    be found first, the match(es) for the second observation will be found second, etc.
    Matching without replacement will generally increase bias.
    Ties are randomly broken when \code{replace==FALSE}
    ---see the \code{ties} option for details.}
  \item{ties}{A logical flag for whether ties should be handled deterministically.  By
    default \code{ties==TRUE}. If, for example, one treated observation
    matches more than one control observation, the matched dataset will
    include the multiple matched control observations and the matched data
    will be weighted to reflect the multiple matches.  The sum of the
    weighted observations will still equal the original number of
    observations. If \code{ties==FALSE}, ties will be randomly broken.
    \emph{If the dataset is large and there are many ties, setting
    \code{ties=FALSE} often results in a large speedup.} Whether two
    potential matches are close enough to be considered tied, is
    controlled by the \code{distance.tolerance}
    option.}
  \item{CommonSupport}{This logical flag implements the usual procedure
    by which observations outside of the common support of a variable
    (usually the propensity score) across treatment and control groups are
    discarded.  The \code{caliper} option is to
    be preferred to this option because \code{CommonSupport}, consistent
    with the literature, only drops \emph{outliers} and leaves
    \emph{inliers} while the caliper option drops both.
    If \code{CommonSupport==TRUE}, common support will be enforced on
    the first variable in the \code{X} matrix.  Note that dropping
    observations generally changes the quantity being estimated. Use of
    this option renders it impossible to use the returned
    objects \code{index.treated} and \code{index.control} to
    reconstruct the matched dataset.  The returned object \code{mdata}
    will, however, still contain the matched dataset.  Seriously, don't
    use this option; use the \code{caliper} option instead.}
  \item{Weight}{ A scalar for the type of weighting scheme the matching
    algorithm should use when weighting each of the covariates in
    \code{X}.  The default value of 1 denotes that weights are equal to
    the inverse of the variances. 2 denotes the Mahalanobis distance
    metric, and 3 denotes that the user will supply a weight matrix
    (\code{Weight.matrix}).  Note that if the user supplies a
    \code{Weight.matrix}, \code{Weight} will be automatically set to be
    equal to 3.}
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
  \item{weights}{ A vector the same length as \code{Y} which
    provides observation specific weights.}
  \item{Var.calc}{ A scalar for the variance estimate
    that should be used.  By default \code{Var.calc=0} which means that
    homoscedasticity is assumed.  For values of  \code{Var.calc > 0},
    robust variances are calculated using \code{Var.calc} matches. }
  \item{sample}{ A logical flag for whether the population or sample
    variance is returned. }
  \item{distance.tolerance}{This is a scalar which is used to determine
    if distances between two observations are different from zero.  Values
    less than \code{distance.tolerance} are deemed to be equal to zero.
    This option can be used to perform a type of optimal matching}
  \item{tolerance}{ This is a scalar which is used to determine
    numerical tolerances.  This option is used by numerical routines
    such as those used to determine if a matrix is singular.}
  \item{restrict}{A matrix which restricts the possible matches.  This
    matrix has one row for each restriction and three
    columns.  The first two columns contain the two observation numbers
    which are to be restricted (for example 4 and 20), and the third
    column is the restriction imposed on the observation-pair.
    Negative numbers in the third column imply that the two observations
    cannot be matched under any circumstances, and positive numbers are
    passed on as the distance between the two observations for the
    matching algorithm.  The most commonly used positive restriction is
    \code{0} which implies that the two observations will always
    be matched.  \cr

    Exclusion restrictions are even more common.  For example, if we want
    to exclude the observation pair 4 and 20 and
    the pair 6 and 55 from being matched, the restrict matrix would be:
    \code{restrict=rbind(c(4,20,-1),c(6,55,-1))}}
  \item{match.out}{ The return object from a previous call to
    \code{Match}.  If this object is provided, then \code{Match} will
    use the matches found by the previous invocation of the function.
    Hence, \code{Match} will run faster.  This is
    useful when the treatment does not vary across calls to
    \code{Match} and one wants to use the same set of matches as found
    before.  This often occurs when one is trying to estimate the causal
    effect of the same treatment (\code{Tr}) on different outcomes
    (\code{Y}). When using this option, be careful to use the same
    arguments as used for the previous invocation of \code{Match} unless
    you know exactly what you are doing.}
  \item{version}{The version of the code to be used.  The "fast" C/C++
    version of the code does not calculate Abadie-Imbens standard errors.
    Additional speed can be obtained by setting \code{ties=FALSE} or
    \code{replace=FALSE} if the dataset is large and/or has many ties.
    The "legacy" version of the code does not make a call to an optimized
    C/C++ library and is included only for historical compatibility.  The
    "fast" version of the code is significantly faster than the "standard"
    version for large datasets, and the "legacy" version is much slower
    than either of the other two.}
}
\details{

  This function is intended to be used in conjunction with the
  \code{MatchBalance} function which checks if the results of this
  function have actually achieved balance.  The results of this function
  can be summarized by a call to the \code{\link{summary.Match}}
  function. If one wants to do propensity score matching, one should estimate the
  propensity model before calling \code{Match}, and then place the
  fitted values in the \code{X} matrix---see the provided example.  \cr

  The \code{\link{GenMatch}} function can be used to \emph{automatically
  find balance} by the use of a genetic search algorithm which determines
  the optimal weight to give each covariate. The object returned by
  \code{\link{GenMatch}} can be supplied to the \code{Weight.matrix}
  option of \code{Match} to obtain estimates.\cr

  \code{Match} is often much faster with large datasets if
  \code{ties=FALSE} or \code{replace=FALSE}---i.e., if matching is done
  by randomly breaking ties or without replacement.  Also see the
  \code{\link{Matchby}} function.  It provides a wrapper for
  \code{Match} which is much faster for large datasets when it can be
  used.\cr

  Three demos are included: \code{GerberGreenImai}, \code{DehejiaWahba},
  and \code{AbadieImbens}.  These can be run by calling the
  \code{\link{demo}} function such as by \code{demo(DehejiaWahba)}. \cr
}

\value{
  \item{est }{The estimated average causal effect.}
  \item{se }{The Abadie-Imbens standard error.  This standard error has
    correct coverage if \code{X} consists of either covariates or a known
    propensity score because it takes into account the uncertainty of the
    matching procedure.  If an estimated propensity score is used, the
    uncertainty involved in its estimation is not accounted for although
    the uncertainty of the matching procedure itself still is.}
  \item{est.noadj }{The estimated average causal effect without any
    \code{BiasAdjust}.  If \code{BiasAdjust} is not requested, this is the
    same as \code{est}.}
  \item{se.standard }{The usual standard error.  This is the standard error
    calculated on the matched data using the usual method of calculating
    the difference of means (between treated and control) weighted by the
    observation weights provided by \code{weights}.  Note that the
    standard error provided by \code{se} takes into account the uncertainty
    of the matching procedure while \code{se.standard} does not.  Neither
    \code{se} nor \code{se.standard} take into account the uncertainty of
    estimating a propensity score.  \code{se.standard} does
    not take into account any \code{BiasAdjust}.  Summary of both types
    of standard error results can be requested by setting the
    \code{full=TRUE} flag when using the \code{\link{summary.Match}}
    function on the object returned by \code{Match}.}
  \item{se.cond }{The conditional standard error.  The practitioner
    should not generally use this.}
  \item{mdata }{A list which contains the matched datasets produced by
    \code{Match}.  Three datasets are included in this list: \code{Y},
    \code{Tr} and \code{X}.}
  \item{index.treated }{A vector containing the observation numbers from
    the original dataset for the treated observations in the
    matched dataset.  This index in conjunction with \code{index.control}
    can be used to recover the matched dataset produced by
    \code{Match}.  For example, the \code{X} matrix used by \code{Match}
    can be recovered by
    \code{rbind(X[index.treated,],X[index.control,])}. The user should
    generally just examine the output of \code{mdata}.}
  \item{index.control }{A vector containing the observation numbers from
    the original data for the control observations in the
    matched data.  This index in conjunction with \code{index.treated}
    can be used to recover the matched dataset produced by
    \code{Match}.  For example, the \code{X} matrix used by \code{Match}
    can be recovered by
    \code{rbind(X[index.treated,],X[index.control,])}. The user should
    generally just examine the output of \code{mdata}.}
  \item{index.dropped}{A vector containing the observation numbers from
    the original data which were dropped (if any) in the matched dataset
    because of various options such as \code{caliper} and
    \code{exact}.  If no observations were dropped, this
    index will be \code{NULL}.}
  \item{weights}{A vector of weights.  There is one weight for each
    matched-pair in the matched dataset.  If all of the observations had a
    weight of 1 on input, then each matched-pair will have a weight of 1
    on output if there are no ties.}
  \item{orig.nobs }{The original number of observations in the dataset.}
  \item{orig.wnobs }{The original number of weighted observations in the
    dataset.}
  \item{orig.treated.nobs}{The original number of treated observations (unweighted).}
  \item{nobs }{The number of observations in the matched dataset.}
  \item{wnobs }{The number of weighted observations in the matched dataset.}
  \item{caliper }{The \code{caliper} which was used.}
  \item{ecaliper }{The size of the enforced caliper on the scale of the
    \code{X} variables.  This object has the same length as the number of
    covariates in \code{X}.}
  \item{exact}{The value of the \code{exact} function argument.}
  \item{ndrops}{The number of weighted observations which were dropped
    either because of caliper or exact matching.  This number, unlike
    \code{ndrops.matches}, takes into account observation specific
    weights which the user may have provided via the \code{weights}
    argument.}
  \item{ndrops.matches}{The number of matches which were dropped either because of
    caliper or exact matching.}
}
\references{
  Sekhon, Jasjeet S. 2011.  "Multivariate and Propensity Score
  Matching Software with Automated Balance Optimization.''
  \emph{Journal of Statistical Software} 42(7): 1-52.
  \doi{10.18637/jss.v042.i07}

  Diamond, Alexis and Jasjeet S. Sekhon. 2013. "Genetic
  Matching for Estimating Causal Effects: A General Multivariate
  Matching Method for Achieving Balance in Observational Studies.''
  \emph{Review of Economics and Statistics}. 95 (3): 932--945.
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
\seealso{ Also see \code{\link{summary.Match}},
  \code{\link{GenMatch}},
  \code{\link{MatchBalance}},
  \code{\link{Matchby}},
  \code{\link{balanceUV}},
  \code{\link{qqstats}}, \code{\link{ks.boot}},
  \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
}
\examples{

# Replication of Dehejia and Wahba psid3 model
#
# Dehejia, Rajeev and Sadek Wahba. 1999.``Causal Effects in
# Non-Experimental Studies: Re-Evaluating the Evaluation of Training
# Programs.''Journal of the American Statistical Association 94 (448):
# 1053-1062.

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
# Estimating the treatment effect on the treated (the "estimand" option defaults to ATT).
#
rr  <- Match(Y=Y, Tr=Tr, X=X, M=1);
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
%  LocalWords:  nodegr rr nboots nmc mb Neyman Tr CommonSupport Var var Matchby
%  LocalWords:  est orig nobs Guido ealso qqstats educ nonparametric rep es ks
%  LocalWords:  inliers Mahalanobis
