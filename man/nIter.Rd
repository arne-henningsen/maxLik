\name{nIter}
\alias{nIter}
\alias{nIter.default}
\title{Return number of iterations for iterative models} 
\description{
  Returns the number of iterations for iterative models.  The default
  method assumes presence of a component \code{$iterations} in {x}.
}
\usage{
nIter(x, \dots)
\method{nIter}{default}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, or a result of maximisation, such as
    created by \code{\link{maxLik}} or \code{\link{maxNR}}}
  \item{\dots}{further arguments for methods}
}
\details{
  This is a generic function.  The default method returns the component
  \code{x$iterations}.
}
\value{
  numeric, number of iterations
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{maxLik}}, \code{\link{maxNR}}
  }
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Estimate with numeric gradient and numeric Hessian
a <- maxNR(loglik, start=1)
nIter(a)
}
\keyword{methods}
