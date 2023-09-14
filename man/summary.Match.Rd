\name{summary.Match}
\alias{summary.Match}
\alias{print.summary.Match}
\title{Summarizing output from Match}
\description{
  \code{\link{summary}} method for class \code{\link{Match}}
}
\usage{
  \method{summary}{Match}(object, ... , full=FALSE, digits=5)
}
\arguments{
  \item{object}{An object of class "\code{Match}", usually, a
    result of a call to \code{\link{Match}}.}
  \item{full}{A flag for whether the unadjusted estimates and naive
    standard errors should also be summarized.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{https://www.jsekhon.com}.  
}
\seealso{ Also see \code{\link{Match}}, \code{\link{GenMatch}},
  \code{\link{MatchBalance}}, 
  \code{\link{balanceUV}}, \code{\link{qqstats}}, \code{\link{ks.boot}},
  \code{\link{GerberGreenImai}}, \code{\link{lalonde}}
}
\keyword{htest}

%  LocalWords:  UC url seealso GenMatch MatchBalance balanceUV htest
%  LocalWords:  GerberGreenImai lalonde
