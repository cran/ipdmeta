\name{beta.indices}
\alias{beta.indices}

\title{Shared Components for Mixed-Level Meta-Analysis Functions}

\description{
Gives position for shared coefficients of individual-level and
study-level models of mixed-level meta-analyses.
}

\usage{
beta.indices(ipd.f,ad.f)
}

\arguments{

  \item{ipd.f}{Formula of ipd model as given in \code{mlma}}
  \item{ad.f}{Study-level formula}
}

\value{

List with components:

\tabular{ll}{
 \code{ipd} \tab Positions for full model beta that are in IPD model
 \cr
\code{ad} \tab Positions for full model beta that are in AD model
 \cr
\code{shared.index} \tab Componets of AD coefficients that are shared
with IPD (logical) \cr
\code{labels} \tab Coefficient labels for full beta 
}

}

\seealso{
\code{\link{mlma}}
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

ipd.f <- Surv(time,event)~trt*I(x-x.bar)+trt*x.bar
ad.f <- surv~log(time)+trt*x.bar

beta.indices(ipd.f,ad.f)

}
