\name{coxmcem}
\alias{coxmcem}

\title{MCEM with importance sampling for Cox Mixed Effects Model}

\description{
Estimation for the parameters of a Cox mixed effects model using an EM algorithm with importance sampling to perform the E step.
}

\usage{

coxmcem(

	fixed,
	random,
	data,
	df,
	n.groups,
	max.iter=10,
	min.sample=100,
	mc.step=2,
	est.delta=5/100,
	echo = TRUE,
	init.coef,
	init.vcov
)

}

\arguments{

  \item{fixed}{Formula for fixed effects model as would be given in \code{coxph}}
  \item{random}{Formula for random effects as would be given in \code{coxme}, with explicit removal of intercept term if not to be included, i.e. ~-1}
  \item{data}{Data frame of observed data}
  \item{df}{Degrees of freedom for multivariate T proposal}
  \item{n.groups}{Number of clusters}
  \item{max.iter}{Maximum number of EM iterations}
  \item{min.sample}{Starting sample size for frailties at expectation step}
  \item{mc.step}{Integer which increases MC samples if error is still too large}
  \item{est.delta}{Converge criterion for change in estimates between EM iterations}
  \item{echo}{Logical, whether to print convergence details at each EM iteration}
  \item{init.coef}{Initial values for fixed effect parameters}
  \item{init.vcov}{Intiial variance-covariance matrix parameters}
}


\value{

Returns a list with :

	\item{coef}{Fixed effects estimates}
	\item{cluster}{Importance-weighted mean cluster effects}
	\item{vcov}{Estimated cluster effect variance}
	\item{var}{List with \code{coef} and \code{vcov} which give the variance for the fixed effects and random effect variances (diagonal components only), respectively}
	\item{convergence}{Vector of maximum relative difference of all model parameters between consecutive iterations}
	\item{mc.samples}{Vector of number of draws from joint proposal for cluster effects which are the basis for the importance-weighted averages at the E step}
	\item{sd.loglik}{Standard deviation of the target log-likelihood for the sampled frailties. Useful for monitoring variability in the joint frailties and proposal distribution diffuseness.}
	\item{loglik}{Conditional log-likelihood of full model}
	\item{loglik.fixed}{Log-likelihood of model with no random effects}
	\item{max.weight}{Vector of the iteration's maximum importance weight}
	\item{iterations}{The number of EM iterations performed}
}

\details{

If \code{init.coef} and \code{init.vcov} are not specified, starting values for the regression effects and frailty variance are set to the values in a \code{coxme} or \code{phmm} is used to set initial values, with \code{coxme} tried first. 

The proposal density used for the cluster effects at each E step is a multivariate T distribution. Independent samples are taken with location and scale set to the \code{mean} and \code{Sigma} of \code{coxmcem.initialize}. The \code{mean} is determined from the frailties of a PPL-based estimation implemented by \code{coxme}. If this fails then the MCEM method of \cite{phmm} is used. The value for \code{scale} is a diagonal matrix with entries that are the reciprocal of the effective sample size for each frailty which is the total number of observations divided by the total number of frailties in the model.


Weights are based on the ratio of the penalized partial likelihood and the proposal density times an upper bound. The upper bound is the maximum of the PPL to the proposal among the E-step samples. These weights are then normalized so as to sum to 1.

Convergence is declared reached if the maximum relative change has been less than \code{est.delta} for three consecutive iterations.

The sample size changes with the result of the CV check of the relative change in the estimates, comparing the most recent three consecuritve differences to the previous three consecutive differences. If the current iteration CV is greater than the last the sample size increases by \eqn{N+\frac{N}{mc.step}}{N + N/mc.step}.

}

\seealso{
\code{\link{coxme}}, \code{\link{coxph}}, \code{\link{phmm}}
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

set.seed(3232)

data(cancer.rats)

fit.rats <- 

coxmcem(Surv(time,event)~rx,~(1|litter),n.groups=50,data=cancer.rats,
		max.iter=10,
		min.sample=300,
		mc.step=2,
		est.delta=.01,
		df=50
		)

}
