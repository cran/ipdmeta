\name{converged}
\alias{converged}
\alias{converged,ipdlme-method}
\title{Extract the convergence status for an \code{ipdlme} object}
\description{
  Extract the convergence status for the fit of an IPD LME model.
}
\usage{
converged(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
  Logical value indicating whether convergence was achieved.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

converged(ipdlme(n, Y, S2))
}

\keyword{models}
