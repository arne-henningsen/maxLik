\name{vcov.maxLik}
\alias{vcov.maxLik}

\title{Variance Covariance Matrix of maxLik objects}
\description{
   Extract variance-covariance matrices of objects of class \code{\link{maxLik}}.
}
\usage{
   \method{vcov}{maxLik}( object, eigentol=1e-12, ... )
}
\arguments{
   \item{object}{an object of class \code{probit} or \code{maxLik}.}
   \item{eigentol}{
     nonzero print limit on the range of the absolute values of the
     hessian.  Specifically, define:

     absEig <- eigen(hessian(object), symmetric=TRUE)[['values']]

     Then compute and print t values, p values, etc. only if min(absEig)
     > (eigentol * max(absEig)).  
   }
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen,
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
vcov(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
vcov(a)
}

\keyword{methods}
