\name{maxControl}
\alias{maxControl}
\title{
  Optimization control options
}
\description{
  This is the preferred way to supply optimization control parameters.
  It takes the parameters, performs consistency checks, and returns the
  control structure.
}
\usage{
maxControl(...)
}
\arguments{
  \item{\dots}{
    list of optimization parameters.  The supported parameters and their
    default values are
    \itemize{
      \item{tol}{1e-8, stopping condition 
	for \code{\link{maxNR}} and related optimizers.
	Stop if the absolute difference
	between successive iterations is less than \code{tol}, return
	\code{code=2}.}
      \item{reltol}{sqrt(.Machine$double.eps), relative convergence
	tolerance
	for \code{\link{maxNR}} and related optimizers.
	The algorithm stops if
        it is unable to increase the value by a factor of 'reltol *
        (abs(val) + reltol)' at a step.}
      \item{gradtol}{1e-6, stopping condition
	for \code{\link{maxNR}} and related optimizers.
	Stop if norm of the gradient is
	less than \code{gradtol}, return \code{code=1}.}
      \item{steptol}{1e-10, stopping/error condition
	for \code{\link{maxNR}} and related optimizers.
	If \code{qac == "stephalving"} and the quadratic
	approximation leads to a worse instead of a better value, or to
	\code{NA}, the step length
	is halved and a new attempt is made.  If necessary, this procedure is repeated
	until step < \code{steptol}, thereafter \code{code=3} is returned.}
      %
      \item{lambdatol}{1e-6, (for \code{\link{maxNR}} related
	optimizers)
	controls whether Hessian is treated as negative
	definite.  If the
	largest of the eigenvalues of the Hessian is larger than
	\code{-lambdatol}, a suitable diagonal matrix is subtracted from the
	Hessian (quadratic hill-climbing) in order to enforce nagetive
	definiteness.}
      %
      \item{qac}{"stephalving", character, Qadratic Approximation
	Correction for \code{\link{maxNR}} related optimizers.  When the new
	guess is worse than the initial one, program attemts to correct it:
	\code{"stephalving"} decreases the
	step but keeps the direction,
	\code{"marquardt"} uses
	Marquardt (1963) method by decreasing the step length while also
	moving closer to the pure gradient direction.  It may be faster and
	more robust choice in areas where quadratic approximation behaves poorly.}
      \item{qrtol}{1e-10, QR-decomposition tolerance
	for \code{\link{maxNR}} related optimizers.
      }
      \item{lambda0}{1e-2, a positive numeric, initial correction term
	for \cite{Marquardt (1963)} correction in
	\code{\link{maxNR}}-related optimizers}
      \item{lambdaStep}{2, how much the \cite{Marquardt (1963)}
	correction is decreased/increased at
	successful/unsuccesful step
	for \code{\link{maxNR}} related optimizers}
      \item{maxLambda}{1e12, maximum allowed correction term
	for \code{\link{maxNR}} related optimizers.
	If exceeded, the
	algorithm exits with return code 3.}
      %
      \item{iterlim}{150L, stopping condition.  Stop if more than \code{iterlim}
	iterations performed.  Note that \sQuote{iteration} may mean
	different things for different optimzers.}
      \item{printLevel}{0L, the level of verbosity.  Larger values print
	more information.  Depends on the optimizer.}
    }
}
}
\value{
  object of class 'MaxControl', the control structure to supply for the
  optimization routines.
}
\author{
  Ott Toomet \email{otoomet@gmail.com}
}
\note{
  Several of the control parameters can also be supplied directly to the
  optimization routines.
}
\examples{
## Optimize quadratic form t(D) %*% W %*% D with p.d. weight matrix,
## s.t. constraints sum(D) = 1
quadForm <- function(D) {
   return(-t(D) \%*\% W \%*\% D)
}
eps <- 0.1
W <- diag(3) + matrix(runif(9), 3, 3)*eps
D <- rep(1/3, 3)
                        # initial values
library(maxLik)
res <- maxNR(quadForm, start=D,
             control=maxControl(print.level=2,
                     qac="marquardt", lambda0=1)
             )
print(summary(res))
}
\keyword{utilities}
