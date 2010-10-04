\name{ipd.meta.power}
\alias{ipd.meta.power}

\title{IPD meta-analysis covariate-treatment power estimate}

\description{
Estimates the individual-patient meta-analysis power of a Wald test of
interaction in a mixed effects model for a continuous or binary
outcome with only study-level summaries (no patient-level data).
}

\usage{
ipd.meta.power(
y0,
y1,
var0,
var1,
x0,
x1,
s20,
s21,
n0,
n1,
interaction,
alpha=.05,
logOR=TRUE,
level=95)
}

\arguments{

All vector entries are grouped by study.

  \item{y0}{vector of control group response}
  \item{y1}{vector of treatment  group response}
  \item{var0}{vector of variances of patient-level responses for control group}
  \item{var1}{vector of variances of patient-level responses for treatment group}
  \item{x0}{vector of covariate means for control group}
  \item{x1}{vector of covariate means for control group}
  \item{s20}{vector of covariate sample variances for control group}
  \item{s21}{vector of covariate sample variances for treatment group}
  \item{n0}{vector of number of control subjects}
  \item{n1}{vector of number of treated subjects}

 \item{interaction}{scalar, supposed alternative interaction effect}
 \item{alpha}{scalar significance level of Wald test (two-sided)}
 \item{logOR}{logical, whether the effect size is logit odds ratio}
 \item{level}{0 to 100 level for confidence interval of variance/power estimate}
}

\details{

The estimator applies to an IPD generalized linear mixed model with a
logit or idenity link, corresponding to a binary or continuous
outcome, respectively. The supposed IPD GLMM has an intervention,
covariate and covariate-treatment interaction term. Random effects by
study for baseline and treatment effect are assumed. 

The study-level data for a continuous response has effect sizes that
are the mean difference between groups;  summary responses are the treatment arm means. For a binary outcome, the effect size is the log-odds ratio and the response for each treatment arm is the log odds of the mean proportion. When the effect size is the log-odds the residual variance \code{sigma} is always one as is the case for a GLM model with logit link and no overdispersion.

When logOR is \code{FALSE} it is assumed that the patient-level response is a continuous variable. 

}

\value{

A list with the following named components:

\tabular{ll}{
 \code{esimated.power} \tab The estmated IPD meta-analysis interactive effect power \cr
 \code{power.lower} \tab Lower bound for \code{level} CI \cr
 \code{power.upper} \tab Upper bound for \code{level} CI \cr
 \code{estimated.se} \tab Estimated standard error of IPD meta-analysis interaction effect \cr
 \code{se.lower} \tab Lower bound for \code{level} CI \cr
 \code{se.upper} \tab Upper bound for \code{level} CI \cr
 \code{sigma} \tab Estimate of residual variance \cr
\code{sigma0} \tab Estimate of intercept random effect variance from simple RE meta-analysis with DL estimator \cr
\code{sigma1} \tab Estimate of treatment random effect variance simple RE meta-analysis with DL estimator \cr
\code{level} \tab 0 to 100 level of uncertainty intervals
 }
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

data(poynard)

#CREATE META SUMMARY OF LOG-OR EFFECT SIZE FOR 2-YR BLEEDING EVENTS

binary.mean <- function(x){
n <- length(x)
if(sum(x)==0){
 return(.5/(n+.5))
 }
else if(sum(x)==n){
 return(n/(n+.5))
 }
else{
return(mean(x))
 }
}

logit <- function(x){log(x/(1-x))}

meta <- data.frame(p0=sapply((poynard$bleed0/poynard$n0),binary.mean))

meta$p1 <- sapply((poynard$bleed1/poynard$n1),binary.mean)
meta$y0 <- logit(meta$p0)
meta$y1 <- logit(meta$p1)
meta$var0 <- 1/(meta$p0*(1-meta$p0))
meta$var1 <- 1/(meta$p1*(1-meta$p1))
meta$n0 <- poynard$n0
meta$n1 <- poynard$n1

#ALTERNATIVE HYPOTHESIS FOR AGE-TREATMENT EFFECT
#WITH 10 YEARS CHANGE TO OR TREATMENT EFFECT exp(beta*10)
#EFFECT MODIFIER CHANGES TREATMENT EFFECT BY 30%

beta = log(1.3)/10

ipd.power <- 

ipd.meta.power(
 int=beta,
 y0=meta$y0,y1=meta$y1,
 var0=meta$var0,var1=meta$var1,
 n0=meta$n0,n1=meta$n1,
 x0=poynard$age0,x1=poynard$age1,
 s20=poynard$age.s20,s21=poynard$age.s21
)

ipd.power

}
