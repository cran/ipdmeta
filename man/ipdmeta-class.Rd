\name{ipdlme-class}
\title{Class for representing ipdlme objects}
\docType{class}
% Classes
\alias{ipdlme-class}
% Methods
\alias{coef,ipdlme-method}
\alias{summary,ipdlme-method}
%\alias{confint,ipdlme-method}
\alias{print,ipdlme-method}
\alias{show,ipdlme-method}
\alias{summary,ipdlme-method}
\alias{vcov,ipdlme-method}


\description{The \code{ipdlme} class represents the fit of a linear mixed-effects model of individual patient data meta-analysis of multiple parallel group clinical trials, based on aggregate data estimation methods.
}

\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("ipdlme", ...)}  or via the function
  \code{\link{ipdlme}}.
}

\section{Slots}{
  The class \code{"ipdlme"} represents a linear mixed model to assess effect modification of multiple clinical trials and contains the slots:

  \describe{

    \item{\code{fixef}:}{The vector of the population effect estimates.}

    \item{\code{ranef}:}{The matrix of the study-specific intercept and treatment random effects.}

    \item{\code{vcov.fixef}:}{The variance-covariance matrix for inference with the population effects.}

    \item{\code{vcov.ranef}:}{The variance-covariance matrix for inference with the study random effects.}

    \item{\code{sigma2}:}{Estimate of the residual variance.}

    \item{\code{VarCorr}:}{Covariance-variance matrix for bivariate normal random effects.}

    \item{\code{convergence.trace}:}{A list of the values of each of the effects and variance components at each iteration in the maximization algorith.}

     \item{\code{converged}:}{Logical value indicating whether the convergence criterion was met.}

    \item{\code{n.iter}:}{The total number of iterations used in the optimization algorithm.}

    \item{\code{max.iter}:}{The maximum number of iterations specified for the optimization algorithm.}

    \item{\code{tol}:}{The tolerance level on which the termination rule for convergence is based.}

    \item{\code{df}:}{The degrees of freedom of the model fit.}

  }% summary slots
}

\section{Methods}{

The following methods are extractors for the component that shares the method's name:

  \describe{
\item{\code{fixef}}{\code{signature(object = "ipdlme")}{}}
\item{\code{ranef}}{\code{signature(object = "ipdlme")}{}}
\item{\code{coef}}{\code{signature(object = "ipdlme")}{}}
\item{\code{vcov}}{\code{signature(object = "ipdlme")}{}}
\item{\code{Var}}{\code{signature(object = "ipdlme")}{}}
\item{\code{sigma2}}{\code{signature(object = "ipdlme")}{}}
\item{\code{vcov.fixef}}{\code{signature(object = "ipdlme")}{}}
\item{\code{vcov.ranef}}{\code{signature(object = "ipdlme")}{}}
\item{\code{convergence}}{\code{signature(object = "ipdlme")}{}}
\item{\code{converged}}{\code{signature(object = "ipdlme")}{}}
\item{\code{n.iter}}{\code{signature(object = "ipdlme")}{}}
\item{\code{tol}}{\code{signature(object = "ipdlme")}{}}
\item{\code{max.iter}}{\code{signature(object = "ipdlme")}{}}
	}

\describe{				      
    \item{\code{print}}{\code{signature(x = "ipdlme")}: print information about
      the fitted model.  }
    \item{\code{show}}{\code{signature(object = "ipdlme")}: Same as the
      \code{print} method.}
      \item{\code{confint}}{\code{signature(object = "ipdlme",parm, level = 0.95, ...)}Returns the specified confidence interval for all the population parameters.}
      \item{\code{plot}}{\code{signature(x = "ipdlme",y,...)}: Displays a forest plot of the study intercept and treatment effects with the option of user-defined labels for the studies.}
\item{\code{summary}}{\code{signature(object = "ipdlme")}:Summary table of standard error and Wald tests for the population effects. A list of the study random effects and estimates of the variance components are also displayed.}
 }

}
\seealso{\code{\link{ipdlme}}}

\examples{
data(regress_chol)

metafit <- ipdlme(n,Y,S2)

converged(metafit)

summary(metafit)

confint(metafit)


}

\keyword{classes}
