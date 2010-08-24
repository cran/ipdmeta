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
   init.beta0,
   init.beta,
   error=1/1000
)
}

\arguments{

  \item{ipd.formula}{Formula for patient-level Cox model as would be supplied to \code{coxph}}
  \item{meta.formula}{Formula for study-level data with only the covariates that are not included in patient-level model, i.e. \code{surv~log(time)}}
  \item{ipd.data}{Patient-level data frame}
  \item{meta.data}{Study-level data frame containing survival estimates by study and treatment group}
  \item{sigma2}{Estimates variances for aggregate survival estimates of \code{meta.data}}
  \item{study.group.interaction}{Factor that is the study and treatment group interaction for \code{meta.data}}
  \item{init.beta0}{Starting values for \code{beta0} additional regression parameters for the aggregate data according to the GLM survival model.}
  \item{init.beta}{Starting values for \code{beta} the shared parameters of the patient-level and study-level.}
  \item{error}{Newton-Raphson stopping rule criterion based on change in log-likelihood}
}


\details{

The estimation is based on likelihood methods for the combined likelihood contributions of the individual and aggregate data. The individual-level model is a Cox proportional hazards model with fixed effects of covariates. The aggregate data are assumed to be a set of survival estimates and their variance (e.g. KM estimates) for each cluster of the meta dataset, i.e. treatment group within study. 

The aggregate model is a multivariate GLM on the log-negative-log transform of the survival estimates. The mean of the transformed survival is linear in the covariate model of \code{meta.formula}. The transformed estimates are assumed to have an MVN distribution where the covariance-variance matrix is regarded as known and is determined by the variances \code{sigma2} and the correlation between surival estimates within the same cluster.

The formula argument for \code{ipd.formula} specifies the set of shared covariates, i.e. all model components save intercept term and factors for the baseline hazard. These covariates share the same variable name in \code{ipd.data} and \code{meta.data}.

The formula argument \code{meta.formula} will have components which model the way that time is related to the survival estimate. This function should be the first term following \code{surv~}. It is assumed that all the components of \code{ipd.formula} are also covariates for the aggregate data so that these dod not nee to be specified in the \code{meta.formula}.

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
    Surv(time,event)~trt,surv~log(time),
    ipd.data,
    meta.data,
    meta.data$sigma2,
    meta.data$sub.group,
    init.beta0=c(1,0),
    init.beta=-.5
    )

fit$coef		#MODEL FIT

sqrt(diag(fit$var))	#STANDARD ERROR

}
