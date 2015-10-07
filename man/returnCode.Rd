\name{returnCode}
\alias{returnCode}
\alias{returnCode.default}
\alias{returnCode.maxLik}
\alias{returnMessage}
\alias{returnMessage.default}
\alias{returnMessage.maxim}
\alias{returnMessage.maxLik}
\title{Success or failure of the optimization}
\description{
  These function extract success or failure information from optimization objects.
  The \code{returnCode} gives a numeric code, and \code{returnMessage} a
  brief description about the success or
  failure of the optimization, and point to the problems occured (see
  documentation for the
  corresponding functions).
}
\usage{
returnCode(x, ...)
\method{returnCode}{default}(x, ...)
\method{returnCode}{maxLik}(x, ...)
returnMessage(x, ...)
\method{returnMessage}{maxim}(x, ...)
\method{returnMessage}{maxLik}(x, ...)
}
\arguments{
  \item{x}{object, usually an optimization result}
  \item{...}{further arguments for other methods}
}
\details{
  \code{returnMessage} and \code{returnCode} are a generic functions, with methods
  for various optimisation algorithms.
  The message should either describe
  the convergence (stopping condition),
  or the problem.  
}
\value{
  Integer for \code{returnCode}, character for \code{returnMessage}.
  Different optimization routines may define it in a different way.
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}}, \code{\link{maxBFGS}}}
\examples{
## maximise the exponential bell
f1 <- function(x) exp(-x^2)
a <- maxNR(f1, start=2)
returnCode(a) # should be success (1 or 2)
returnMessage(a)
## Now try to maximise log() function
a <- maxNR(log, start=2)
returnCode(a) # should give a failure (4)
returnMessage(a)
}
\keyword{methods}
\keyword{utilities}
