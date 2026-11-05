\name{summary.maxLik}
\alias{summary.maxLik}
\alias{coef.summary.maxLik}
\title{summary the Maximum-Likelihood estimation}
\description{
  Summary the Maximum-Likelihood estimation including standard errors
  and t-values.
}
\usage{
\method{summary}{maxLik}(object, eigentol=1e-12, ... )
\method{coef}{summary.maxLik}(object, \dots)
}

\arguments{
  \item{object}{
    object of class 'maxLik', or 'summary.maxLik',
    usually a result from
    Maximum-Likelihood estimation.
  }
  \item{eigentol}{
    The standard errors are only calculated if the ratio of the smallest
    and largest eigenvalue of the Hessian matrix is less than
    \dQuote{eigentol}.  Otherwise the Hessian is treated as singular.
  }
  \item{\ldots}{currently not used.}
}
\value{
  An object of class \code{summary.maxLik} summarizing the estimation
  results.  This information can be extracted with the following
  methods:
  \describe{
    \item{\code{\link[=coef.summary.maxLik]{coef}}}{numeric matrix, the
      first column contains the parameter estimates, the second the
      standard errors, the third the t-values and the fourth the
      corresponding probabilities.  The standard errors are \code{Inf}
      if the Hessian is treated as singular, and the t-values and
      probabilities are \code{NA} for fixed parameters.}
    \item{\code{\link[=logLik.maxLik]{logLik}}}{the log-likelihood value
      in the maximum; its \dQuote{df} attribute gives the number of free
      parameters.}
    \item{\code{\link{maximType}}}{type of maximization.}
    \item{\code{\link{nIter}}}{number of iterations.}
    \item{\code{\link{returnCode}}}{code of success.}
    \item{\code{\link{returnMessage}}}{a short message describing the
      code.}
    \item{\code{\link{activePar}}}{logical vector indicating which
      parameters were free, i.e. not treated as constants.}
  }
  The variance-covariance matrix, the standard errors and the number of
  observations are extracted from the fitted \code{maxLik} object
  itself, using \code{\link[=vcov.maxLik]{vcov}},
  \code{\link[=stdEr.maxLik]{stdEr}} and \code{\link[=nobs.maxLik]{nobs}}.

  Information about a constrained maximization, if any, can currently
  only be extracted directly (with \code{\$}): the returned object
  contains a \code{constraints} component in this case, otherwise
  \code{NULL}.  It is a list passed on from the fitted object,
  describing the type of constrained optimization, the number of outer
  iterations and the barrier value, plus a return code and a message
  when these are available (they are missing e.g. for
  \code{constrOptim}).
}
\author{Ott Toomet, Arne Henningsen}
\seealso{
  \code{\link{maxLik}} for maximum likelihood estimation,
  \code{\link{confint}} for confidence intervals, and \code{\link{tidy}}
  and \code{\link{glance}} for alternative quick summaries of the ML
  results. 
}
\examples{
## ML estimation of exponential distribution:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, control=list(printLevel=2))
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1, control=list(printLevel=2))
summary(a)
}
\keyword{models}
