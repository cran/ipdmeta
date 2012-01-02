\name{tol}
\alias{tol}
\alias{tol,ipdlme-method}
\title{Extract the convergence criterion for an \code{ipdlme} object}
\description{
  Extract the convergence criterion for the fit of an IPD LME model.
}
\usage{
tol(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
Numeric convergence threshold for max percentage change in population and study effect estimates that must be met to declare convergence.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

tol(ipdlme(n, Y, S2))
}

\keyword{models}
