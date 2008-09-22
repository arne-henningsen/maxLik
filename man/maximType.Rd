\name{maximType}
\alias{maximType}
\alias{maximType.default}
\alias{maximType.maxim}
\alias{maximType.MLEstimate}
\title{Type of Minimization/Maximization}
\description{
  Returns the type of optimisation.  It should be returned
  by the optimisation routine.
}
\usage{
  maximType(x)
}
\arguments{
  \item{x}{object of class 'maxim' or another object which
    involves numerical optimisation.
  }
}
\value{
  A text message, describing the involved optimisation algorithm
}
\author{Ott Toomet, \email{otoomet@ut.ee}}
\seealso{\code{\link{maxNR}}}
\examples{
## maximise two-dimensional exponential hat.  Maximum is at c(2,1):
f <- function(a) exp(-(a[1] - 2)^2 - (a[2] - 1)^2)
m <- maxNR(f, start=c(0,0))
summary(m)
maximType(m)
## Now use BFGS maximisation.
m <- maxBFGS(f, start=c(0,0))
summary(m)
maximType(m)
}
\keyword{optimize}
\keyword{methods}
