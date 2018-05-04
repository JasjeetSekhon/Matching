\name{summary.balanceUV}
\alias{summary.balanceUV}
\title{Summarizing output from balanceUV}
\description{
  \code{\link{summary}} method for class \code{\link{balanceUV}}
}
\usage{
  \method{summary}{balanceUV}(object, ..., digits=5)
}
\arguments{
  \item{object}{An object of class "\code{balanceUV}", usually, a
    result of a call to \code{\link{balanceUV}}.}
  \item{digits}{The number of significant digits that should be
    displayed.}
  \item{...}{Other options for the generic summary function.}  
}
\author{
  Jasjeet S. Sekhon, UC Berkeley, \email{sekhon@berkeley.edu},
  \url{http://sekhon.berkeley.edu/}.  
}
\seealso{
  Also see \code{\link{balanceUV}}, 
  \code{\link{Match}}, \code{\link{GenMatch}},
  \code{\link{MatchBalance}}, \code{\link{qqstats}}, \code{\link{ks.boot}},
  \code{\link{GerberGreenImai}},
  \code{\link{lalonde}}
}
\keyword{htest}

%  LocalWords:  balanceUV UC url seealso GenMatch MatchBalance htest
%  LocalWords:  GerberGreenImai lalonde
