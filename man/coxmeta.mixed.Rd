\name{coxmeta.mixed}
\alias{coxmeta.mixed}

\title{Mixed data mixed effects meta-analysis with survival data}

\description{
Estimates fixed and random effects for general frailty structure in a mixed data meta-analysis of censored data.
}


\usage{
coxmeta.mixed(

   ipd.formula,
   meta.formula,
   random.formula,
   ipd.data,
   meta.data,
   ipd.groups,
   meta.groups,
   sigma2,
   study.group.interaction,
   p,
   max.iter,
   min.sample,
   est.delta,
   mc.step,
   df,
   echo = TRUE,
   max.sample,
   converge.index
)
}

\arguments{

  \item{ipd.formula}{Formula for patient-level Cox model as would be supplied to \code{coxph}}
  \item{meta.formula}{Formula for study-level data with only the covariates that are not included in patient-level model, i.e. \code{surv~log(time)}}
  \item{random.formula}{Formula of random effects conditional on group factor as would be specified in \code{coxme}, i.e. \code{~(1+trt|study)}}
  \item{ipd.data}{Patient-level data frame}
  \item{meta.data}{Study-level data frame containing survival estimates by study and treatment group}
  \item{ipd.groups}{Number of clusters in patient-level dataset}
  \item{meta.groups}{Number of clusters in study-level dataset}
  \item{sigma2}{Variances for aggregate survival estimates of \code{meta.data}}
  \item{study.group.interaction}{Factor that is the study and treatment group interaction for \code{meta.data}}
  \item{p}{Number of shared fixed effect parameters, i.e. covariates in \code{ipd.formula}}
  \item{max.iter}{Maximum iterations of MCEM algorithm}
  \item{min.sample}{Starting sample size for samples at E step} 
    \item{est.delta}{Stopping criterion threshold based on standardized relative differences for model parameters}
  \item{mc.step}{Determines how sample size is increased to reduce MC error}
  \item{df}{Degrees of freedom for proposal multivariate T distribution}
  \item{echo}{Logical, if TRUE prints maximum weight and standardized difference between iteration estimates of model parameters}
  \item{max.sample}{If given, then a fixed escalation is used so that if max.iter is reached, final sample size is max.sample}
  \item{converge.index}{Index of which parameters to include in the convergence assessment. Default is all parameters}
}


\value{

List of 

\item{iterations}{Number of MCEM iterations}
\item{est.converge}{Maximal difference in parameter estimates between consecutive iterations relative to their standard error}
\item{max.weight}{Maximum importance weight at E step}
\item{sd.loglik}{Standard deviaiton of each importance-sample conditional Log-likelihood, which can be helpful for monitoring the performance of the proposal distribution}
\item{loglik}{Importance-sample average Conditional Log-likelihood at each iteration}
\item{coef}{Final estimated fixed effects}
\item{vcov}{Final estimated frailty variances}
\item{cluster}{Mean frailty effects at final MCEM iteration}
\item{var}{List of variance of model parameters with components \code{coef} and \code{vcov}}

}

\details{

An MCEM algorithm is used to obtain the model estimates. The expectation step utilizes importance sampling based on a multivariate T distribution for the joint frailties. 

\code{meta.data} are the survival estimates from a KM plot by treatment group. The data.frame should also contain the time measure corresponding to each estimate and the same covariates as in \code{ipd.formula} but defined at the aggregate-level. For patient-level covariates, the treatment-within-study sample average is used.

The stopping rule for the EM algorithm is when the maximum of the last three standardized maximum differences are below \code{est.delta} or when \code{max.iter} has been reached.

The MC error is monitored by the CV of the maximal standardized parameter differences, based on the three most recent iterations. If the CV of the current iteration is greater than the last, the sample size is increased according to \eqn{N+\frac{N}{mc.step}}{N+N/mc.step}.

The stopping rule for the M step Newton-Rapshon procedure is based on the relative difference of consecutive log-likelihood, stopping when this difference is <.0001.

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{
\code{\link{coxph}}, \code{\link{coxme}}, \code{\link{coxmeta.fixed}}
}

\examples{

data(ipd.data)
data(meta.data)
set.seed(401)

#UNIVARIATE FRAILTY BY STUDY FACTOR

fit <- coxmeta.mixed(
    Surv(time,event)~trt,surv~log(time),~(1|group),
    ipd.data,
    meta.data,
    ipd.groups=8,meta.groups=2,
    meta.data$sigma2,
    meta.data$sub.group,
    p=1,max.iter=5,
    min.sample=200,est.delta=.05,mc.step=1.5,df=10
    )

fit$coef				#MODEL FIT
sqrt(diag(fit$var$coef))		#STANDARD ERROR


sqrt(diag(fit$vcov))	#ESTIMATED FRAILTY STANDARD DEVIATION


}
