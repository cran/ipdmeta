\name{mlma}
\alias{mlma}

\title{Mixed-Level Meta-Analysis for Survival Data}

\description{
Mixed-level meta-analysis of a time-to-event endpoint. Combines
evidence from individual-level studies and from study-level trials
where aggregate evidence is in the form of survival proportions along
a survival curve, e.g. surv = c(.9,.5,.2) at months = c(1,6,12).

Fixed effects that share the same name in the study-level and
individual-level model formula are treated as having a common
effect. This allows flexibility in separating study-level and
patient-level effects by using distinct term labels.

MLMA allows for a general random effects structure; multivariate
normal frailties.
}


\usage{
mlma(
   ipd.formula,
   meta.formula,
   random.formula,
   ipd.data,
   meta.data,
   ipd.groups,
   meta.groups,
   sigma2,
   study.group.interaction,
   max.iter=5,
   min.sample=300,
   est.delta=.1,
   mc.step=2,
   df=30,
   echo = TRUE,
   max.sample = NULL,
   converge.index,
   error=1/100,
   fixed=FALSE
)
}

\arguments{

  \item{ipd.formula}{Formula for patient-level Cox model as would be
      supplied to \code{coxph}}
  \item{meta.formula}{Formula for study-level data
      \code{surv~log(time)+trt*x}}
  \item{random.formula}{Formula of random effects conditional on group factor as would be specified in \code{coxme}, i.e. \code{~(1+trt|study)}}
  \item{ipd.data}{Patient-level data frame}
  \item{meta.data}{Study-level data frame containing survival estimates by study and treatment group}
  \item{ipd.groups}{Number of clusters in patient-level dataset}
  \item{meta.groups}{Number of clusters in study-level dataset}
  \item{sigma2}{Variances for aggregate survival estimates of \code{meta.data}}
  \item{study.group.interaction}{Factor that is the study and treatment group interaction for \code{meta.data}}
  \item{max.iter}{Maximum iterations of MCEM algorithm}
  \item{min.sample}{Starting sample size for samples at E step} 
    \item{est.delta}{Stopping criterion threshold based on standardized relative differences for model parameters}
  \item{mc.step}{Determines how sample size is increased to reduce MC error}
  \item{df}{Degrees of freedom for proposal multivariate T distribution}
  \item{echo}{Logical, if TRUE prints maximum weight and standardized difference between iteration estimates of model parameters}
  \item{max.sample}{If given, then a fixed escalation is used so that if max.iter is reached, final sample size is max.sample}
  \item{converge.index}{Index of which parameters to include in the
      convergence assessment. Default is all parameters}
  \item{error}{MLE stopping rule; MLE considered reached when change in log-likelihood for NR maximization is
      less than \code{error}}
  \item{fixed}{Logical, should fixed effects model be used}
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

The estimation is based on likelihood methods for the combined
likelihood contributions of the individual- and trial-level data. The
individual-level model is a Cox proportional hazards mixed model (PHMM). The
study-level model is a GLMM on survival proportions (KM estimates) for
each treatment group within study. The link is the
log-complementary-log transform. When all outcomes across studies and
types follow the PHMM, the log-negative-log survival is linear in the
covariates of the individual-level model.

Conditional on the frailty terms, the transformed survival estimates
are assumed to have an MVN distribution where the covariance-variance
matrix is regarded as known and is determined by the KM survival
estimate variances \code{sigma2} and the correlation between survival estimates within the same cluster.

Shared parameters between the individual-level and study-level models
are determined by common term.labels for \code{ipd.formula} and
\code{meta.formula}. So, for example, if \code{ipd.formula} was
\code{Surv(time,event)~trt*x} and \code{meta.formula} was
\code{surv~log(time)+trt*x} then the shared terms are
\code{trt},\code{x} and \code{trt:x} and these would be combined. If
instead \code{meta.formula} was \code{surv~log(time)+trt*x.bar}, then
only \code{trt} would be common, and only study level data would be
used to estimate \code{x.bar} and \code{trt:x.bar} factors.

When random effects are included, an MCEM algorithm is used to obtain
the full model estimates. The expectation step utilizes importance sampling based on a multivariate T distribution for the joint frailties. 

The stopping rule for the EM algorithm is when the maximum of the last three standardized maximum differences are below \code{est.delta} or when \code{max.iter} has been reached.

The MC error is monitored by the CV of the maximal standardized parameter differences, based on the three most recent iterations. If the CV of the current iteration is greater than the last, the sample size is increased according to \eqn{N+\frac{N}{mc.step}}{N+N/mc.step}.

The stopping rule for the M step Newton-Rapshon procedure is based on the relative difference of consecutive log-likelihood, stopping when this difference is <\code{error}

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{
\code{\link{coxph}}, \code{\link{coxme}}
}

\examples{

#MIXED-LEVEL META-ANALYSIS SURVIVAL DATA

data(ipd.data)
example(ipd.data)

data(meta.data)
example(meta.data)

set.seed(401)

#FIXED EFFECTS

fit <- mlma(
    Surv(time,event)~trt,surv~log(time)+trt,
    ipd.data=ipd.data,
    meta.data=meta.data,
    ipd.groups=8,meta.groups=2,
    study.group=meta.data$sub.group,
    sigma2=meta.data$sigma2,
    fixed=TRUE
    )

fit$coef				#MODEL FIT
sqrt(diag(fit$var))	          	#STANDARD ERROR


#MIXED EFFECTS, BASELINE FRAILTY BY GROUP

fit <- mlma(
    Surv(time,event)~trt,surv~log(time)+trt,~(1|group),
    ipd.data,
    meta.data,
    ipd.groups=8,meta.groups=2,
    study.group=meta.data$sub.group,
    sigma2=meta.data$sigma2,
    max.iter=3
    )

fit$coef				#MODEL FIT
sqrt(diag(fit$var$coef))		#STANDARD ERROR

#ECOLOGICAL BIAS
#MIXED-LEVEL WITH ECOLOGICAL BIAS
#MIXED EFFECTS MLMA SEPARATING STUDY-LEVEL AND INDIVIDUAL-LEVEL X EFFECTS
#SIMULATION EFFECTS TO GENERATE MIXED.BIASED

beta = array(c(0,-.5,-.2,0,-.4,0))
names(beta) <- c("int","trt","x","x.bar","trt and x","trt and x.bar")

data(mixed.biased)

fit <- mlma(
     Surv(time,event)~trt*I(x-x.bar)+trt*x.bar,surv~log(time)+trt*x.bar,
     ~(1|group),
    mixed.biased$ipd,
    mixed.biased$meta,
    ipd.groups=5,meta.groups=5,
    mixed.biased$meta$sigma2,
    mixed.biased$meta$sub.group,
    max.iter=3
    )

fit$coef				#MODEL FIT
sqrt(diag(fit$var$coef))		#STANDARD ERROR


sqrt(diag(fit$vcov))	#ESTIMATED FRAILTY STANDARD DEVIATION


#95\% CI FOR EQUIVALENCE OF STUDY- AND INDIVIDUAL-LEVEL X (MAIN EFFECT)
C <- c(0,-1,1,0,0,0,0)
ci(C,fit$coef,fit$var$coef)

#95\% FOR EQUIVALENCE OF STUDY- AND INDIVIDUAL-LEVEL X-TRT (INTERACTION)
C <- c(0,0,0,-1,1,0,0)
ci(C,fit$coef,fit$var$coef)

}
