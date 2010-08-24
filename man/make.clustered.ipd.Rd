\name{make.clustered.ipd}
\alias{make.clustered.ipd}

\title{Generates a Mixed Effects IPD meta dataset of Weibull event times}

\description{
Generates a random ipd datset for a covariate-treatment Weibull IPD mixed effects meta-analysis. 
}

\usage{

make.clustered.ipd(n,beta,D,censor.rate,tau2)

}

\arguments{

  \item{n}{vector of group sizes}
  \item{beta}{Fixed effects parameters (intercept, trt, x, interaction)}
  \item{D}{Matrix of bivariate frailty covariance-variance}
  \item{censor.rate}{Non-informative target censoring rate}
  \item{tau2}{Between-study variance for normal covariate}
}

\details{

The covariate x is normal with marginal distribution N(0,1+tau2).

The treatment indicator is Bernoulli(1/2).

The \code{beta} parameters are fixed and the study random effects are drawn from a MVN with varcov = D. If these components are zero than the corresponding study effect is set to zero. The shape of the Weibull is drawn from a Uniform(.5,3) and is the same for all studies.

Censoring is a non-informative exponential random variable. The rate is set to be the rate at which the probability of the censoring time being less than the median of the data set is \code{censor.rate}.

}

\value{

List of \code{ipd}, \code{D}, \code{shape} and \code{censor.rate}.

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

make.clustered.ipd(c(5,10,20),rep(0,4),D=diag(c(.1,.2)),.1,1/2)

}
