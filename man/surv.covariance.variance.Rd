\name{surv.covariance.variance}
\alias{surv.covariance.variance}

\title{Construct covariance-variance matrix for set of survival estimates}

\description{
Constructs the covariance-variance matrix for the transformed survival estimates of a meta-analysis (multiple studies, treatment groups and survival estimates). The GLM framework has the link function g(x) = log(-log(x)).
}

\usage{
surv.covariance.variance(surv,sigma2,group)
}

\arguments{

  \item{surv}{vector of survival estimates ordered by group}
  \item{sigma2}{vector of variances for surv estimates (on survival scale)}
  \item{group}{vector of grouping factor which should be treatment and study interaction}

}


\value{

Block-diagonal covariance-variance matrix.
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{\code{\link{bdiag}}}
