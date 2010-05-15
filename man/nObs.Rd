\name{nObs}
\alias{nObs}
\alias{nObs.default}
\alias{nObs.lm}
\alias{nObs.maxLik}
\title{Return number of observations for statistical models} 
\description{
  Returns number of observations for statistical models.  The default
  method assumes presence of a component \code{param$nObs} in \code{x}.
}
\usage{
nObs(x, \dots)
\method{nObs}{default}(x, \dots)
\method{nObs}{lm}(x, \dots)
\method{nObs}{maxLik}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, such as created by \code{\link{lm}}}
  \item{\dots}{further arguments for methods}
}
\details{
  This is a generic function.  The default method returns the component
  \code{x$param$nObs}.  The \code{lm}-method is based on
  qr-decomposition, in the same way as the does \code{\link{summary.lm}}.

   The \code{nObs} method for objects of class \code{"maxLik"}
   can return the number of observations only
   if \code{\link{maxLik}} was called
   with argument \code{grad} equal to a gradient function or
   (if no gradient function is specified)
   argument \code{logLik} equal to a log-likelihood function
   that return the gradients or log-likelihood values, respectively,
   for each observation.
}
\value{
  numeric, number of observations
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{nParam}}}
\examples{
# Construct a simple OLS regression:
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
nObs(m)

## fit a normal distribution by ML
# generate a variable from normally distributed random numbers
x <- rnorm( 100, 1, 2 )
# log likelihood function (for individual observations)
llf <- function( param ) {
   return( dnorm( x, mean = param[ 1 ], sd = param[ 2 ], log = TRUE ) )
}
## ML method
ml <- maxLik( llf, start = c( mu = 0, sigma = 1 ) )
# return number of onservations
nObs( ml )
}
\keyword{methods}
