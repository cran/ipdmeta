\name{MaxIter}
\alias{MaxIter}
\alias{MaxIter,ipdlme-method}
\title{Extract the maximum number of iterations specified for an \code{ipdlme} object}
\description{
  Extract the maximum number of iterations specified for the fit of an IPD LME model.
}
\usage{
MaxIter(object,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Setting for the maximum possible iterations performed in fitting the IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

MaxIter(ipdlme(n, Y, S2))
}

\keyword{models}
