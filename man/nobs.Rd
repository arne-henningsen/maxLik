\name{nObs.maxLik}
\alias{nObs.maxLik}
\alias{nobs,maxLik-method}
\title{Number of Observations}
\description{
  Returns the number of observations for statistical models,
  estimated by Maximum Likelihood using \code{\link{maxLik}}.
  \code{nobs} users the standard generic from \pkg{stats4}, \code{nObs}
  (deprecated)
  uses the S3 generic from \pkg{miscTools}.
}
\usage{
\method{nObs}{maxLik}(x, \dots)
\S4method{nobs}{maxLik}(object, \dots)
}
\arguments{
  \item{object}{a statistical model estimated by Maximum Likelihood
      using \code{\link{maxLik}}.}
  \item{x}{a statistical model estimated by Maximum Likelihood
      using \code{\link{maxLik}}.}
  \item{\dots}{further arguments (currently ignored).}
}
\details{
   The \code{nobs} method for \dQuote{maxLik} objects
   can return the number of observations only if log-likelihood function
   (or the gradient) returns values by individual observation.
}

\value{
  numeric, number of observations
}

\author{Arne Henningsen, Ott Toomet}

\seealso{
  \code{\link[stats]{nobs}},
  \code{\link[miscTools]{nObs}}, \code{\link{maxLik}},
  \code{\link{nParam}}.
}

\examples{
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
nobs( ml )
}
\keyword{methods}
