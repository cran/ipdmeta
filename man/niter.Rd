\name{n.iter}
\alias{n.iter}
\alias{n.iter,ipdlme-method}
\title{Extract the number of iterations for an \code{ipdlme} object}
\description{
  Extract the number of iterations for the fit of an IPD LME model.
}
\usage{
n.iter(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
 Number of iterations of maximization alogrithm perfomed in fitting IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

n.iter(ipdlme(n, Y, S2))
}

\keyword{models}
