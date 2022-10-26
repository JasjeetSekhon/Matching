\name{GerberGreenImai}
\alias{GerberGreenImai}
\docType{data}
\title{Gerber and Green Dataset used by Imai}
\description{

  This is the dataset used by Imai (2005) to replicate and evaluate the
  field experiment done by Gerber and Green (2000). The accompanying
  demo replicates Imai's propensity score model which is then used to
  estimate the causal effect of get-out-the-vote telephone calls on
  turnout. 
  
}
\usage{data(GerberGreenImai)}
\format{
  A data frame with 10829 observations on the following 26 variables.
  \describe{
    \item{PERSONS }{Number persons in household}
    \item{WARD    }{Ward of residence}
    \item{QUESTION}{Asked to commit to voting}
    \item{MAILGRP }{Sent mail}
    \item{PHONEGRP}{Phone batch #1}
    \item{PERSNGRP}{Personal contact attempted}
    \item{APPEAL  }{Content of message}
    \item{CONTACT }{Personal contact occurred}
    \item{MAILINGS}{Number of mailings sent}
    \item{AGE     }{Age of respondent}
    \item{MAJORPTY}{Democratic or Republican}
    \item{VOTE96.0}{Abstained in 1996}
    \item{VOTE96.1}{Voted in 1996}
    \item{MAILCALL}{Phone batch #2}
    \item{VOTED98 }{Voted in 1998}
    \item{PHNSCRPT}{Script read to phone respondents}
    \item{DIS.MC  }{Contacted by phone in batch #2}
    \item{DIS.PHN }{Contacted by phone in batch #1}
    \item{PHN.C   }{Contacted by phone}
    \item{PHNTRT1 }{Phone contact attempted (no blood or blood/civic)}
    \item{PHNTRT2 }{Phone contact attempted (no blood)}
    \item{PHN.C1  }{Contact occurred in phntrt1}
    \item{PHN.C2  }{Contact occurred in phntrt2}
    \item{NEW     }{New voter}
    \item{phone   }{Contacted by phone}
    \item{AGE2    }{Age squared}
  }
}
\details{

  The demo provided, entitled \code{GerberGreenImai}, uses Imai's
  propensity score model to estimate the causal effect of
  get-out-the-vote telephone calls on turnout. The propensity score
  model fails to balance age.

}
\references{

  Gerber, Alan S. and Donald P. Green. 2000. ``The Effects of
  Canvassing, Telephone Calls, and Direct Mail on Voter Turnout: A Field
  Experiment.'' \emph{American Political Science Review} 94: 653-663.

  Gerber, Alan S. and Donald P. Green. 2005. ``Correction to Gerber
  and Green (2000), replication of disputed findings, and reply to
  Imai (2005).''  \emph{American Political Science Review} 99:
  301-313.

  Imai, Kosuke. 2005. ``Do Get-Out-The-Vote Calls Reduce Turnout?
  The Importance of Statistical Methods for Field Experiments.'' 
  \emph{American Political Science Review} 99: 283-300.

  Hansen, Ben B. Hansen and Jake Bowers. forthcoming. ``Attributing
  Effects to a Cluster Randomized Get-Out-The-Vote Campaign.''
  \emph{Journal of the American Statistical Association}.
} 
\seealso{ Also see \code{\link{Match}} and \code{\link{MatchBalance}},
  \code{\link{GenMatch}}, \code{\link{balanceUV}}, \code{\link{ks.boot}}
  \code{\link{lalonde}}
}
\keyword{datasets}


%  LocalWords:  GerberGreenImai docType Dataset dataset Imai's MAILGRP PHONEGRP
%  LocalWords:  PERSNGRP MAJORPTY MAILCALL PHNSCRPT PHN PHNTRT phntrt emph
%  LocalWords:  Kosuke seealso MatchBalance GenMatch balanceUV
%  LocalWords:  lalonde datasets
