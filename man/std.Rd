\name{std}
\alias{std}
\alias{std.default}
\alias{std.lm}

\title{Standard deviations}
\description{
   Extract standard deviations from estimated models.
}
\usage{
  std(x, ...)
\method{std}{default}(x, \dots)
\method{std}{lm}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, such as created by \code{\link{lm}}}
  \item{\dots}{further arguments for methods}
}
\details{
  \code{std} is a generic function with methods for objects of "lm"
  class.  The default method returns the square root of the diagonal of
  the variance-covariance matrix.
}
\value{
  numeric, the estimated standard errors of the coefficients.  
}
\author{
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link{vcov}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
std(a)
}

\keyword{methods}
