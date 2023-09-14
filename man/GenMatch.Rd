\name{GenMatch}
\alias{GenMatch}
\title{Genetic Matching}
\description{
  This function finds optimal balance using multivariate matching where
  a genetic search algorithm determines the weight each covariate is
  given.  Balance is determined by examining cumulative probability
  distribution functions of a variety of standardized statistics.  By
  default, these statistics include t-tests and Kolmogorov-Smirnov
  tests. A variety of descriptive statistics based on empirical-QQ
  (eQQ) plots can also be used or any user provided measure of balance.
  The statistics are not used to conduct formal hypothesis tests,
  because no measure of balance is a monotonic function of bias and
  because balance should be maximized without limit. The object
  returned by \code{GenMatch} can be supplied to the \code{\link{Match}}
  function (via the \code{Weight.matrix} option) to obtain causal
  estimates.  \code{GenMatch} uses \code{\link[rgenoud]{genoud}} to
  perform the genetic search.  Using the \code{cluster} option, one may
  use multiple computers, CPUs or cores to perform parallel
  computations.
}
\usage{
GenMatch(Tr, X, BalanceMatrix=X, estimand="ATT", M=1, weights=NULL,
         pop.size = 100, max.generations=100,
         wait.generations=4, hard.generation.limit=FALSE,
         starting.values=rep(1,ncol(X)),
         fit.func="pvals",
         MemoryMatrix=TRUE,
         exact=NULL, caliper=NULL, replace=TRUE, ties=TRUE,
         CommonSupport=FALSE, nboots=0, ks=TRUE, verbose=FALSE,
         distance.tolerance=1e-05,
         tolerance=sqrt(.Machine$double.eps),
         min.weight=0, max.weight=1000,
         Domains=NULL, print.level=2,
         project.path=NULL,
         paired=TRUE, loss=1,
         data.type.integer=FALSE,
         restrict=NULL,
         cluster=FALSE, balance=TRUE, ...)
}
\arguments{
  \item{Tr}{ A vector indicating the observations which are in the
    treatment regime and those which are not.  This can either be a
    logical vector or a real vector where 0 denotes control and 1 denotes
    treatment.}
  \item{X}{ A matrix containing the variables we wish to match on.
    This matrix may contain the actual observed covariates or the
    propensity score or a combination of both.}
  \item{BalanceMatrix}{ A matrix containing the variables we wish
    to achieve balance on.  This is by default equal to \code{X}, but it can
    in principle be a matrix which contains more or less variables than
    \code{X} or variables which are transformed in various ways.  See
    the examples.}
  \item{estimand}{ A character string for the estimand.  The default
    estimand is "ATT", the sample average treatment effect for the
    treated. "ATE" is the sample average treatment effect, and
    "ATC" is the sample average treatment effect for the controls.}
  \item{M}{A scalar for the number of matches which should be
    found. The default is one-to-one matching. Also see the \code{ties}
    option.}
  \item{weights}{ A vector the same length as \code{Y} which
    provides observation specific weights.}
  \item{pop.size}{Population Size.  This is the number of individuals
    \code{\link[rgenoud]{genoud}} uses to solve the optimization problem.
    The theorems proving that genetic algorithms find good solutions are
    asymptotic in population size.  Therefore, it is important that this value not
    be small.  See \code{\link[rgenoud]{genoud}} for more details.}
  \item{max.generations}{ Maximum Generations.  This is the maximum
    number of generations that \code{\link[rgenoud]{genoud}} will run when
    optimizing.  This is a \emph{soft} limit.  The maximum generation
    limit will be binding only if \code{hard.generation.limit} has been
    set equal to \emph{TRUE}.  Otherwise, \code{wait.generations} controls
    when optimization stops. See \code{\link[rgenoud]{genoud}} for more
    details.}
  \item{wait.generations}{If there is no improvement in the objective
    function in this number of generations, optimization will stop.  The
    other options controlling termination are \code{max.generations} and
    \code{hard.generation.limit}.}
  \item{hard.generation.limit}{ This logical variable determines if the
    \code{max.generations} variable is a binding constraint.  If
    \code{hard.generation.limit} is \emph{FALSE}, then
    the algorithm may exceed the \code{max.generations}
    count if the objective function has improved within a given number of
    generations (determined by \code{wait.generations}).}
  \item{starting.values}{
    This vector's length is equal to the number of variables in \code{X}.  This
    vector contains the starting weights each of the variables is
    given. The \code{starting.values} vector is a way for the user
    to insert \emph{one} individual into the starting population.
    \code{\link[rgenoud]{genoud}} will randomly create the other individuals.  These values
    correspond to the diagonal of the \code{Weight.matrix} as described
    in detail in the \code{\link{Match}} function.}
  \item{fit.func}{The balance metric \code{GenMatch} should optimize.
    The user may choose from the following or provide a function:\cr
    \code{pvals}: maximize the p.values from (paired) t-tests and
    Kolmogorov-Smirnov tests conducted for each column in
    \code{BalanceMatrix}.  Lexical optimization is conducted---see the
    \code{loss} option for details.\cr
    \code{qqmean.mean}: calculate the mean standardized difference in the eQQ
    plot for each variable.  Minimize the mean of these differences across
    variables.\cr
    \code{qqmean.max}:  calculate the mean standardized difference in the eQQ
    plot for each variable.  Minimize the maximum of these differences across
    variables.  Lexical optimization is conducted.\cr
    \code{qqmedian.mean}: calculate the median standardized difference in the eQQ
    plot for each variable.  Minimize the median of these differences across
    variables.\cr
    \code{qqmedian.max}:  calculate the median standardized difference in the eQQ
    plot for each variable.  Minimize the maximum of these differences across
    variables.  Lexical optimization is conducted.\cr
    \code{qqmax.mean}: calculate the maximum standardized difference in the eQQ
    plot for each variable.  Minimize the mean of these differences across
    variables.\cr
    \code{qqmax.max}:  calculate the maximum standardized difference in the eQQ
    plot for each variable.  Minimize the maximum of these differences across
    variables.  Lexical optimization is conducted.\cr
    Users may provide their own \code{fit.func}. The name of the user
    provided function should not be backquoted or quoted.  This function needs
    to return a fit value that will be minimized, by lexical
    optimization if more than one fit value is returned.  The function
    should expect two arguments.  The first being the \code{matches} object
    returned by \code{GenMatch}---see
    below.  And the second being a matrix which contains the variables to
    be balanced---i.e., the \code{BalanceMatrix} the user provided to
    \code{GenMatch}. For an example see
    \url{https://www.jsekhon.com}.}
  \item{MemoryMatrix}{
    This variable controls if \code{\link[rgenoud]{genoud}} sets up a memory matrix.  Such a
    matrix ensures that \code{\link[rgenoud]{genoud}} will request the fitness evaluation
    of a given set of parameters only once. The variable may be
    \emph{TRUE} or \emph{FALSE}.  If it is \emph{FALSE}, \code{\link[rgenoud]{genoud}}
    will be aggressive in
    conserving memory.  The most significant negative implication of
    this variable being set to \emph{FALSE} is that \code{\link[rgenoud]{genoud}} will no
    longer maintain a memory
    matrix of all evaluated individuals.  Therefore, \code{\link[rgenoud]{genoud}} may request
    evaluations which it has  previously requested.  When
    the number variables in \code{X} is large, the memory matrix
    consumes a large amount of RAM.\cr
    
    \code{\link[rgenoud]{genoud}}'s memory matrix will require \emph{significantly} less
    memory if the user sets \code{hard.generation.limit} equal
    to \emph{TRUE}.  Doing this is a good way of conserving
    memory while still making use of the memory matrix structure.}
  \item{exact}{ A logical scalar or vector for whether exact matching
    should be done.  If a logical scalar is
    provided, that logical value is applied to all covariates in
    \code{X}.  If a logical vector is provided, a logical value should
    be provided for each covariate in \code{X}. Using a logical vector
    allows the user to specify exact matching for some but not other
    variables.  When exact matches are not found, observations are
    dropped.  \code{distance.tolerance} determines what is considered to
    be an exact match. The \code{exact} option takes precedence over the
    \code{caliper} option.  Obviously, if \code{exact} matching is done
    using \emph{all} of the covariates, one should not be using
    \code{GenMatch} unless the \code{distance.tolerance} has been set
    unusually high.}
  \item{caliper}{ A scalar or vector denoting the caliper(s) which
    should be used when matching.  A caliper is the distance which is
    acceptable for any match.  Observations which are outside of the
    caliper are dropped. If a scalar caliper is provided, this caliper is
    used for all covariates in \code{X}.  If a vector of calipers is
    provided, a caliper value should be provided for each covariate in
    \code{X}. The caliper is interpreted to be in standardized units.  For
    example, \code{caliper=.25} means that all matches not equal to or
    within .25 standard deviations of each covariate in \code{X} are
    dropped.  The \code{ecaliper} object which is returned by
    \code{GenMatch} shows the enforced caliper on the scale of the
    \code{X} variables. Note that dropping observations generally changes
    the quantity being estimated.}
  \item{replace}{A logical flag for whether matching should be done with
    replacement.  Note that if \code{FALSE}, the order of matches
    generally matters.  Matches will be found in the same order as the
    data are sorted.  Thus, the match(es) for the first observation will
    be found first, the match(es) for the second observation will be found second, etc.
    Matching without replacement will generally increase bias.
    Ties are randomly broken when \code{replace==FALSE}---see the
    \code{ties} option for details.}
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
    observations generally changes the quantity being estimated.  Use of
    this option renders it impossible to use the returned
    object \code{matches} to reconstruct the matched dataset.
    Seriously, don't use this option; use the \code{caliper} option instead.}
  \item{nboots}{The number of bootstrap samples to be run for the
    \code{ks} test.  By default this option is set to zero so no
    bootstraps are done.  See \code{\link{ks.boot}} for additional
    details.}
  \item{ks}{ A logical flag for if the univariate bootstrap
    Kolmogorov-Smirnov (KS) test should be calculated.  If the ks option
    is set to true, the univariate KS test is calculated for all
    non-dichotomous variables.  The bootstrap KS test is consistent even
    for non-continuous variables.  By default, the bootstrap KS test is
    not used. To change this see the \code{nboots} option. If a given
    variable is dichotomous, a t-test is used even if the KS test is requested.  See
    \code{\link{ks.boot}} for additional details.}
  \item{verbose}{A logical flag for whether details of each
    fitness evaluation should be printed.  Verbose is set to FALSE if
    the \code{cluster} option is used.}
  \item{distance.tolerance}{This is a scalar which is used to determine
    if distances between two observations are different from zero.  Values
    less than \code{distance.tolerance} are deemed to be equal to zero.
    This option can be used to perform a type of optimal matching.}
  \item{tolerance}{ This is a scalar which is used to determine
    numerical tolerances.  This option is used by numerical routines
    such as those used to determine if a matrix is singular.}
  \item{min.weight}{This is the minimum weight any variable may be
    given.}  
  \item{max.weight}{This is the maximum weight any variable may be
    given.}
  \item{Domains}{This is a \code{ncol(X)} \eqn{\times 2}{*2} matrix.
    The first column is the lower bound, and the second column is the
    upper bound for each variable over which \code{\link[rgenoud]{genoud}} will
    search for weights.  If the user does not provide this matrix, the
    bounds for each variable will be determined by the \code{min.weight}
    and \code{max.weight} options.}
  \item{print.level}{
    This option controls the level of printing.  There
    are four possible levels: 0 (minimal printing), 1 (normal), 2
    (detailed), and 3 (debug).  If level 2 is selected, \code{GenMatch} will
    print details about the population at each generation, including the
    best individual found so far. If debug
    level printing is requested, details of the \code{\link[rgenoud]{genoud}}
    population are printed in the "genoud.pro" file which is located in
    the temporary \code{R} directory returned by the \code{\link{tempdir}}
    function.  See the \code{project.path} option for more details.
    Because \code{GenMatch} runs may take a long time, it is important for the
    user to receive feedback.  Hence, print level 2 has been set as the
    default.}
  \item{project.path}{ This is the path of the
    \code{\link[rgenoud]{genoud}} project file.  By default no file is
    produced unless \code{print.level=3}.  In that case,
    \code{\link[rgenoud]{genoud}} places its output in a file called
    "genoud.pro" located in the temporary directory provided by
    \code{\link{tempdir}}.  If a file path is provided to the
    \code{project.path} option, a file will be created regardless of the
    \code{print.level}. The behavior of the project file, however, will
    depend on the \code{print.level} chosen.  If the \code{print.level}
    variable is set to 1, then the project file is rewritten after each
    generation.  Therefore, only the currently fully completed generation
    is included in the file.  If the \code{print.level} variable is set to
    2 or higher, then each new generation is simply appended to the
    project file. No project file is generated for
    \code{print.level=0}.}
  \item{paired}{A flag for whether the paired \code{\link{t.test}} should be
    used when determining balance.}
  \item{loss}{The loss function to be optimized.  The default value, \code{1},
    implies "lexical" optimization: all of the balance statistics will
    be sorted from the most discrepant to the least and weights will be
    picked which minimize the maximum discrepancy. If multiple sets of
    weights result in the same maximum discrepancy, then the second
    largest discrepancy is examined to choose the best weights.  The
    processes continues iteratively until ties are broken.  \cr

    If the value of \code{2} is used, then only the maximum discrepancy
    is examined.  This was the default behavior prior to version 1.0.  The
    user may also pass in any function she desires. Note that the
    option 1 corresponds to the \code{\link{sort}} function and option 2
    to the \code{\link{min}} function.  Any user specified function
    should expect a vector of balance statistics ("p-values") and it
    should return either a vector of values (in which case "lexical"
    optimization will be done) or a scalar value (which will be
    maximized). Some possible alternative functions are
    \code{\link{mean}} or \code{\link{median}}.}
  \item{data.type.integer}{
    By default, floating-point weights are considered. If this option is
    set to \code{TRUE}, search will be done over integer weights. Note
    that before version 4.1, the default was to use integer weights.}
  \item{restrict}{A matrix which restricts the possible matches. This
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

    Exclusion restriction are even more common.  For example, if we want
    to exclude the observation pair 4 and 20 and the pair 6 and 55 from
    being matched, the restrict matrix would be:
    \code{restrict=rbind(c(4,20,-1),c(6,55,-1))}} \item{cluster}{This
    can either be an object of the 'cluster' class returned by one of
    the \code{\link[parallel]{makeCluster}} commands in the
    \code{parallel} package or a vector of machine names so that
    \code{GenMatch} can setup the cluster automatically. If it is the
    latter, the vector should look like: \cr
    \code{c("localhost","musil","musil","deckard")}.\cr This vector
    would create a cluster with four nodes: one on the localhost another
    on "deckard" and two on the machine named "musil".  Two nodes on a
    given machine make sense if the machine has two or more chips/cores.
    \code{GenMatch} will setup a SOCK cluster by a call to
    \code{\link[parallel]{makePSOCKcluster}}.  This will require the
    user to type in her password for each node as the cluster is by
    default created via \code{ssh}.  One can add on usernames to the
    machine name if it differs from the current shell: "username@musil".
    Other cluster types, such as PVM and MPI, which do not require
    passwords, can be created by directly calling
    \code{\link[parallel]{makeCluster}}, and then passing the returned
    cluster object to \code{GenMatch}. For an example of how to manually
    setup up a cluster with a direct call to
    \code{\link[parallel]{makeCluster}} see
    \url{https://www.jsekhon.com}.  For
    an example of how to get around a firewall by ssh tunneling see:
    \url{https://www.jsekhon.com}.}
    \item{balance}{This logical flag controls if load balancing is done
    across the cluster.  Load balancing can result in better cluster
    utilization; however, increased communication can reduce
    performance.  This option is best used if each individual call to
    \code{\link[Matching]{Match}} takes at least several minutes to
    calculate or if the nodes in the cluster vary significantly in their
    performance. If cluster==FALSE, this option has no effect.}
    \item{...}{Other options which are passed on to
    \code{\link[rgenoud]{genoud}}.}  } \value{ \item{value}{The fit
    values at the solution.  By default, this is a vector of p-values
    sorted from the smallest to the largest.  There will generally be
    twice as many p-values as there are variables in
    \code{BalanceMatrix}, unless there are dichotomous variables in this
    matrix.  There is one p-value for each covariate in
    \code{BalanceMatrix} which is the result of a paired t-test and
    another p-value for each non-dichotomous variable in
    \code{BalanceMatrix} which is the result of a Kolmogorov-Smirnov
    test. Recall that these p-values cannot be interpreted as hypothesis
    tests.  They are simply measures of balance.}  \item{par}{A vector
    of the weights given to each variable in \code{X}.}
    \item{Weight.matrix}{A matrix whose diagonal corresponds to the
    weight given to each variable in \code{X}.  This object corresponds
    to the \code{Weight.matrix} in the \code{\link{Match}} function.}
    \item{matches}{A matrix where the first column contains the row
    numbers of the treated observations in the matched dataset. The
    second column contains the row numbers of the control
    observations. And the third column contains the weight that each
    matched pair is given.  These objects may not correspond
    respectively to the \code{index.treated}, \code{index.control} and
    \code{weights} objects which are returned by \code{\link{Match}}
    because they may be ordered in a different way. Therefore, end users
    should use the objects returned by \code{\link{Match}} because those
    are ordered in the way that users expect.}
  \item{ecaliper }{The
    size of the enforced caliper on the scale of the \code{X} variables.
    This object has the same length as the number of covariates in
    \code{X}.}  } \references{ Sekhon, Jasjeet S. 2011.  "Multivariate
    and Propensity Score Matching Software with Automated Balance
    Optimization.''  \emph{Journal of Statistical Software} 42(7): 1-52.
    \doi{10.18637/jss.v042.i07}
  
  Diamond, Alexis and Jasjeet S. Sekhon. 2013. "Genetic
  Matching for Estimating Causal Effects: A General Multivariate
  Matching Method for Achieving Balance in Observational Studies.''
  \emph{Review of Economics and Statistics}.  95 (3): 932--945.
  \url{https://www.jsekhon.com}

  Sekhon, Jasjeet Singh and Walter R. Mebane, Jr. 1998.  "Genetic
  Optimization Using Derivatives: Theory and Application to Nonlinear
  Models.'' \emph{Political Analysis}, 7: 187-210.
  \url{https://www.jsekhon.com} 
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{https://www.jsekhon.com}.  
}
\seealso{ Also see \code{\link{Match}}, \code{\link{summary.Match}},
  \code{\link{MatchBalance}}, \code{\link[rgenoud]{genoud}},
  \code{\link{balanceUV}}, \code{\link{qqstats}},
  \code{\link{ks.boot}}, \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
} 
\examples{
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
#For details see https://www.jsekhon.com.
#
genout <- GenMatch(Tr=treat, X=X, BalanceMatrix=BalanceMat, estimand="ATE", M=1,
                   pop.size=16, max.generations=10, wait.generations=1)

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

# For more examples see: https://www.jsekhon.com.
}
\keyword{nonparametric}

%  LocalWords:  GenMatch Kolmogorov multinomial CPUs BalanceMatrix estimand ATT
%  LocalWords:  ncol MemoryMatrix nboots ATC rgenoud genoud emph cr ecaliper UC
%  LocalWords:  eqn tempdir rbind makeCluster localhost musil deckard usernames
%  LocalWords:  makePSOCKcluster PVM MPI url dataset seealso MatchBalance cbind
%  LocalWords:  balanceUV GerberGreenImai lalonde hisp nodegr genout
%  LocalWords:  BalanceMat mout mb mv pvals qqmean QQ qqmedian qqmax Smirnov Tr
%  LocalWords:  func CommonSupport min Abadie Imbens boostrap ssh username educ
%  LocalWords:  ealso qqstats nonparametric eQQ rep ks inliers es
