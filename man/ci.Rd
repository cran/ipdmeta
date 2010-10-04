\name{ci}
\alias{ci}

\title{
Confidence Interval for linear combination of parameter estimates.
}

\description{
Returns point-estimate for transformation of linear combination of parameter estimates.
}

\usage{
ci(C,estimates,variance,alpha=.05,f=function(x){exp(x)})
}

\arguments{

\item{C}{Vector or matrix of linear combinations, if matrix
    linear combinations are by row.}
\item{estimates}{Vector of parameter estimates}
\item{variance}{Matrix of variances for estimates}
\item{alpha}{Significance level for normal-based CI}
\item{f}{Function for the parameter transformation at CI endpoints}

} 

\value{

List of transformed point-estimate and CI.

}

\details{

If \code{x} is a matrix, then the \code{ci} function is applied by
row.

}

\author{
S. A. Kovalchik <skoval@ucla.edu>
}

\examples{

data(cars)
fit <- lm(dist~speed,cars)

#AVERAGE DISTANCE GIVEN SPEED OF 10

ci(matrix(c(1,10),1),coef(fit),vcov(fit),f=function(x){x})

}

