\name{bread.maxLik}
\alias{bread.maxLik}

\title{Bread for Sandwich Estimator}
\description{
   Extracting an estimator for the \sQuote{bread} of the sandwich estimator,
   see \code{\link[sandwich]{bread}}.
}
\usage{
\method{bread}{maxLik}( x, ... )
}
\arguments{
   \item{x}{an object of class \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   Matrix, the inverse of the expectation of the second derivative
   (Hessian matrix) of the log-likelihood function
   with respect to the parameters,
   usually equal to the variance covariance matrix of the parameters
   times the number of observations.
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

\seealso{\code{\link[sandwich]{bread}}, \code{\link{maxLik}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t

## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1 )

# Extract the "bread"
library( sandwich )
bread( a )

all.equal( bread( a ), vcov( a ) * nObs( a ) )
}

\keyword{methods}
