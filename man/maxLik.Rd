\name{maxLik}
\alias{maxLik}
\alias{coef.maxLik}
\alias{print.maxLik}

\title{Maximum likelihood estimation}
\description{
  This is just a wrapper for maximisation routines which return object
  of class "maxLik".  Corresponding methods can correctly handle the
  likelihood-specific properties of the estimate including the fact that
  inverse of negative hessian is the variance-covariance matrix.
}
\usage{
maxLik(logLik, grad = NULL, hess = NULL, start, method,
constraints=NULL, ...)
}
\arguments{
  \item{logLik}{log-likelihood function.  Must have the parameter vector
    as the first argument.  Must return either a single log-likelihood
    value or a numeric vector where each component is log-likelihood
    corresponding to individual observations.}
  \item{grad}{gradient of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return either single gradient
    vector with length equal to the number of parameters, or a matrix
    where each row corresponds to gradient vector of individual
    observations.  If \code{NULL}, numeric gradient will be used.}
  \item{hess}{hessian of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return a square matrix.  If
    \code{NULL}, numeric gradient will be used.}
  \item{start}{numeric vector, initial value of parameters.}
  \item{method}{maximisation method, currently either "Newton-Raphson",
    "BFGS", "BFGS-YC",
    "BHHH", "SANN" or "NM" (for Nelder-Mead).  Lower-case letters and
    shortcuts (as 'nr' for Newton-Raphson) allowed.  If missing, a
    suitable method is selected automathically.}
  \item{constraints}{either \code{NULL} for unconstrained maximization
    or a list, specifying the constraints.  See \code{\link{maxBFGS}}.
  }
  \item{\dots}{further arguments for the maximisation routine.}
}
\details{
  \code{maxLik} "support" constrained optimization in the sense that
  constraints are passed further to the underlying optimization
  routines, and suitable default method is selected.  However, no
  attempt is made to correct the resulting variance-covariance matrix.
  Hence, the inference may be wrong.  A corresponding warning is issued
  by the summary method.
}
\value{
  object of class 'maxLik' which inherits from class 'maxim'.
  Components are identical to those of class 'maxim',
  see \code{\link{maxNR}}.
}
\section{Warning}{The constrained maximum likelihood estimation should
  be considered as experimental.  In particular, the variance-covariance
  matrix is not corrected for constrained parameter space.
  }
\author{Ott Toomet \email{otoomet@ut.ee}, Arne Henningsen}
\seealso{\code{\link{maxNR}}, \code{\link{nlm}} and \code{\link{optim}}
  for different non-linear optimisation routines.}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
print( a )
coef( a )
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
print( a )
coef( a )
##
##
## Next, we give an example with vector argument:  Estimate the mean and
## variance of a random normal sample by maximum likelihood
##
loglik <- function(param) {
  mu <- param[1]
  sigma <- param[2]
  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
  ll
}
x <- rnorm(1000, 1, 2) # use mean=1, stdd=2
N <- length(x)
res <- maxLik(loglik, start=c(0,1)) # use 'wrong' start values
print( res )
coef( res )
}
\keyword{optimize}
