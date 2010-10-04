\name{coxmeta.fixed}
\alias{coxmeta.fixed}

\title{Mixed data fixed effects meta-analysis with survival data}

\description{
Combined estimation of parameters of hazard model for a survival meta-analysis when the meta dataset had both patient-level and study-level evidence (mixed data). The fixed effects model assumes homogeneity between studies.
}


\usage{


coxmeta.fixed(

    ipd.formula,
    meta.formula,
    ipd.data,
    meta.data,
    sigma2,
    study.group.interaction,
    beta.ad,
    beta.ipd,
    error=1/1000
)
}

\arguments{

  \item{ipd.formula}{Formula for patient-level Cox model as would be
      supplied to \code{coxph}, \code{Surv(time,event)~trt*x}}
  \item{meta.formula}{Formula for study-level data \code{surv~log(time)+trt*x}}
  \item{ipd.data}{Patient-level data frame}
  \item{meta.data}{Study-level data frame containing survival estimates by study and treatment group}
  \item{sigma2}{Estimates variances for aggregate survival estimates of \code{meta.data}}
  \item{study.group.interaction}{Factor that is the study and treatment group interaction for \code{meta.data}}
  \item{beta.ipd}{Vector of starting values for the individual-level model.}
  \item{beta.ad}{Vector of starting values for the study-level model.}
  \item{error}{Newton-Raphson stopping rule criterion based on change in log-likelihood}
}


\details{

The estimation is based on likelihood methods for the combined likelihood contributions of the individual and aggregate data. The individual-level model is a Cox proportional hazards model with fixed effects of covariates. The aggregate data are assumed to be a set of survival estimates and their variance (e.g. KM estimates) for each cluster of the meta dataset, i.e. treatment group within study. 

The aggregate model is a multivariate GLM on the log-negative-log transform of the survival estimates. The mean of the transformed survival is linear in the covariate model of \code{meta.formula}. The transformed estimates are assumed to have an MVN distribution where the covariance-variance matrix is regarded as known and is determined by the variances \code{sigma2} and the correlation between surival estimates within the same cluster.

Shared parameters between the individual-level and study-level models
are determined by common term.labels for \code{ipd.formula} and
\code{meta.formula}. So, for example, if \code{ipd.formula} was
\code{Surv(time,event)~trt*x} and \code{meta.formula} was
\code{surv~log(time)+trt*x} then the shared terms are
\code{trt},\code{x} and \code{trt:x} and these would be combined. If
instead \code{meta.formula} was \code{surv~log(time)+trt*x.bar}, then
only \code{trt} would be common, and only study level data would be
used to estimate \code{x.bar} and \code{trt:x.bar} factors.
}

\value{

A list with components \code{coef}, \code{var}, \code{loglik}, and a list \code{monitor} which shows the Newton-Raphson changes in \code{coef} and \code{loglik}.

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{
\code{\link{coxmeta.mixed}}
}

\examples{

data(ipd.data)
data(meta.data)

#SHARED COVARIATE ~TRT+X
#SURV BASELINE HAZARD COMPONENT OF MODEL ~LOG(TIME)

fit <- coxmeta.fixed(
    Surv(time,event)~trt,surv~log(time)+trt,
    ipd.data,
    meta.data,
    meta.data$sigma2,
    meta.data$sub.group,
    beta.ipd = runif(1),
    beta.ad = c(1,0,runif(1)),
    )

fit$coef		#MODEL FIT

sqrt(diag(fit$var))	#STANDARD ERROR

###SEPARATION OF STUDY-LEVEL AND PATIENT-LEVEL COVARIATE

ipd.data$x.bar <-
rep(tapply(ipd.data$x,ipd.data$group,mean),table(ipd.data$group))
ipd.data$x.c <- ipd.data$x-ipd.data$x.bar
meta.data$x.bar <- meta.data$x

fit <- coxmeta.fixed(
    Surv(time,event)~trt+x.c+x.bar,surv~log(time)+trt+x.bar,
    ipd.data,
    meta.data,
    meta.data$sigma2,
    meta.data$sub.group,
    beta.ipd = runif(3),
    beta.ad = c(1,0,runif(2)),
    )

fit$coef		#MODEL FIT

sqrt(diag(fit$var))	#STANDARD ERROR


}
