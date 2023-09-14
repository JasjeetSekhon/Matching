\name{summary.ks.boot}
\alias{summary.ks.boot}
\alias{print.summary.ks.boot}
\title{Summarizing output from ks.boot}
\description{
  \code{\link{summary}} method for class \code{\link{ks.boot}}
}
\usage{
  \method{summary}{ks.boot}(object, ..., digits=5)
}
\arguments{
  \item{object}{An object of class "\code{ks.boot}", usually, a
    result of a call to \code{\link{ks.boot}}.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}  
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{https://www.jsekhon.com}.  
}
\seealso{
  Also see \code{\link{ks.boot}}, \code{\link{balanceUV}}, \code{\link{qqstats}},
  \code{\link{Match}}, \code{\link{GenMatch}},
  \code{\link{MatchBalance}}, \code{\link{GerberGreenImai}},
  \code{\link{lalonde}}
}
\keyword{htest}

%  LocalWords:  UC url seealso balanceUV GenMatch MatchBalance htest
%  LocalWords:  GerberGreenImai lalonde
