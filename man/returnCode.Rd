\name{returnCode}
\alias{returnCode}
\alias{returnCode.default}
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
  corresponding function).
}
\usage{
returnCode(x, ...)
\method{returnCode}{default}(x, ...)
returnMessage(x, ...)
\method{returnMessage}{maxim}(x, ...)
\method{returnMessage}{maxLik}(x, ...)
}
\arguments{
  \item{x}{object, usually an estimator, achieved by optimization}
  \item{...}{further arguments for other methods}
}
\details{
  \code{returnMessage} is a generic function, with methods
  for various optimisation algorithms.
  The default methods returns component \code{returnCode}.
  The message should either describe
  the convergence,
  or the problem.  
}
\value{
  Integer, the success code of optimization procedure.  However,
  different optimization routines may define it in a different way.
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
