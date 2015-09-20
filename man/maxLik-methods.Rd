\name{AIC.maxLik}
\alias{AIC.maxLik}
\alias{coef.maxim}
\alias{coef.maxLik}
\alias{stdEr.maxLik}

\title{Methods for the various standard functions}
\description{
  These are methods for the maxLik related objects.  See also the
  documentation for the
  corresponding generic functions
}
\usage{
\method{AIC}{maxLik}(object, \dots, k=2)
\method{coef}{maxim}(object, \dots)
\method{stdEr}{maxLik}(x, \dots)
}
\arguments{
  \item{object}{a \sQuote{maxLik} object (also a
    \sQuote{maxim} object for \code{coef})}
  \item{k}{numeric, the penalty per parameter to be used; the default
    \sQuote{k = 2} is the classical AIC.}
  \item{x}{a \sQuote{maxLik} object}
  \item{\dots}{other arguments for methods}
}
\examples{
## estimate mean and variance of normal random vector
set.seed( 123 )
x <- rnorm(50, 1, 2 )
## log likelihood function.
## Note: 'param' is a vector
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   return(sum(dnorm(x, mean=mu, sd=sigma)))
}
## Estimate it.  Take standard normal as start values
ml <- maxLik(llf, start = c(mu=0, sigma=1) )

coef(ml)
stdEr(ml)
AIC(ml)
}
\keyword{methods}
