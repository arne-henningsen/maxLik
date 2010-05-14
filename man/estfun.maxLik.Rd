\name{estfun.maxLik}
\alias{estfun.maxLik}

\title{Extract Gradients Evaluated at each Observation}
\description{
   Extract the gradients of the log-likelihood function evaluated
   at each observation (\sQuote{Empirical Estimating Function},
   see \code{\link[sandwich]{estfun}}).
}
\usage{
\method{estfun}{maxLik}( x, ... )
}
\arguments{
   \item{x}{an object of class \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   Matrix of gradients of the log-likelihood function
   at the estimated parameter value
   evaluated at each observation
}
\section{Warnings}{
   The \pkg{sandwich} package must be loaded before this method
   can be used.

   This method works only if \code{\link{maxLik}} was called 
   with argument \code{grad} equal to a gradient function or
   (if no gradient function is specified)
   argument \code{logLik} equal to a log-likelihood function
   that return the gradients or log-likelihood values, respectively,
   for each observation.
}
\author{
   Arne Henningsen
}

\seealso{\code{\link[sandwich]{estfun}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t

## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1 )

# Extract the gradients evaluated at each observation
library( sandwich )
estfun( a )

## Estimate with analytic gradient
gradlik <- function(theta) 1/theta - t
b <- maxLik(loglik, gradlik, start=1)
estfun( b )
all.equal( c( estfun( b ) ), gradlik( coef( b ) ) )
}

\keyword{methods}
