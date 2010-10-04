\name{bas.haz}
\alias{bas.haz}

\title{
Baseline hazard for individual-level and study-level from mixed-level
meta-analysis
}

\description{
Returns \code{survfit} object for patient-level data and cumulative
hazard estimates and CI based on study-level estimates.
}

\usage{
bas.haz(
        formula,
        haz.formula,
        data,
        coef,
        var,
        alpha=.05,
        FUN=exp
)
}

\arguments{

\item{formula}{Surv formula for patient-level model as given in
    \code{mlma}.}
\item{haz.formula}{Hazard model component of \code{meta.formula} of
    \code{mlma}}
\item{data}{Patient-level data frame}
\item{coef}{Vector of coefficients output from \code{mlma}}
\item{var}{Matrix of variances for coef estimates from \code{mlma}}
\item{alpha}{Significance level for normal-based CI for study-level
    cumulative hazard}
\item{FUN}{Function for the parameter transformation at CI endpoints}

} 

\value{

List of \code{ipd.survfit} and \code{meta.bas.haz}.

\code{ipd.survfit} is a \code{survfit} object.

\code{meta.bas.haz} is a data frame with \code{est}, \code{lower},
\code{upper} and \code{time}. The time at which the study-level estimates are
evaluated are the event times for the patient-level data set.

}

\details{

The \code{haz.formula} is the component of the study-level model given
in \code{mlma} that corresponds to the log-cumulative hazard,
i.e. \code{~log(t)}.

}

\seealso{
\code{\link{mlma}}, \code{\link{survfit}}
}

\author{
S. A. Kovalchik <skoval@ucla.edu>
}

\examples{

data(ipd.data)
data(meta.data)


fit <- mlma(
     Surv(time,event)~trt,
     surv~-1+log(time),
     random=~(1|group),
     ipd.groups=8,
     meta.groups=2,
     ipd.data=ipd.data,
     meta.data=meta.data,
     sigma2=meta.data$sigma2,
     study.group=meta.data$sub.group,
     max.iter=5,
     min=100
)

#CUMULATIVE HAZARD SUMMARIES

H <- bas.haz(
             Surv(time,event)~trt,
             ~-1+log(t),
             ipd.data,
             fit$coef,
             fit$var$coef,
             )

#RETURN OBJECTS
#PATIENT-LEVEL
H$ipd.survfit

#STUDY-LEVEL
head(H$meta.bas.haz)


#PATIENT-LEVEL BASELINE HAZARD

plot(H$ipd.survfit,fun="cumhaz",ylab="H(t)",bty="n")

#STUDY-LEVEL BASELINE HAZARD WITH 95\% CI

lines(x=H$meta.bas.haz$time,y=H$meta.bas.haz$lower,type="l",lty=2)
lines(x=H$meta.bas.haz$time,y=H$meta.bas.haz$est,type="l",col="blue")
lines(x=H$meta.bas.haz$time,y=H$meta.bas.haz$upper,type="l",lty=2)

legend("topleft",legend=c("Study-level","Patient-level"),lty=1,col=c("blue","black"))
}

