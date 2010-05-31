\name{nParam.maxim}
\alias{nParam.maxim}
\title{Number of model parameters}
\description{
  This function returns the number of model parameters.
}
\usage{
\method{nParam}{maxim}(x, free=FALSE, \dots)
}

\arguments{
  \item{x}{a model returned by a maximisation method from the \pkg{maxLik}
      package.}
  \item{free}{logical, whether to report only the free parameters or the
    total number of parameters (default)}
  \item{\dots}{other arguments for methods}
}
\details{
  Free parameters are the parameters with no equality
  restrictions.  Some parameters may be restricted (e.g. sum of two
  probabilities may be restricted to equal unity).   In this case the
  total number of parameters may depend on the
  normalisation. 
}
\value{
  Number of parameters in the model
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{nObs}} for number of observations}
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
# return number of parameters
nParam( ml )
}
\keyword{methods}
