\name{Var}
\alias{Var}
\alias{Var,ipdlme-method}
\title{Extract the random effect variance parameters}
\description{
  Extract the estimates of study-level baseline and intercept variance parameters for the bivariate normal random effects of the IPD LME model.
}
\usage{
Var(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
  The matrix of the two-by-two covariance-variance for the random effects.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

Var(ipdlme(n, Y, S2))
}

\keyword{models}
