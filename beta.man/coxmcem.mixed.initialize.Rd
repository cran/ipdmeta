\name{coxmcem.mixed.initialize}
\alias{coxmcem.mixed.initialize}

\title{Initialization for Mixed-Level Meta-Analysis T Frailty Proposal}

\description{
Returns patient-level starting values for the multivariate T proposal.
}

\usage{
coxmcem.mixed.initialize(
		   coxme.formula,
		   phmm.formula,
		   data,
		   init.coef,
		   init.vcov)

}

\arguments{

  \item{coxme.formula}{Formula for \code{coxme}}
  \item{phmm.formula}{Formula for \code{phmm}}
  \item{data}{Patient-level data frame}
  \item{init.coef}{Supplied vector for fixed effects}
  \item{init.vcov}{Supplied matrix for frailty variances}

}


\value{

List of \code{mean}, \code{Sigma}  (location and scale of T) and
\code{beta}, \code{D}, the fixed and variance parameters.

}

\details{

First tries \code{coxme}, but if estimation is not possible reverts to
\code{phmm}.

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{\code{\link{beta.indices}}}
