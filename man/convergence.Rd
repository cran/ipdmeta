\name{convergence}
\alias{convergence}
\alias{convergence,ipdlme-method}
\title{Extract the convergence trace for an \code{ipdlme} object}
\description{
  Extract the convergence trace for the fit of an IPD LME model.
}
\usage{
convergence(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
 List of the trace for every parameters of the IPD LME fitted model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

convergence(ipdlme(n, Y, S2))
}

\keyword{models}
