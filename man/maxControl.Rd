\name{MaxControl-class}
\Rdversion{1.1}
\docType{class}
\alias{MaxControl-class}
\alias{maxControl}
\alias{maxControl,MaxControl-method}
\alias{maxControl,missing-method}
\alias{maxControl,maxim-method}
\alias{show,MaxControl-method}

\title{Class \code{"MaxControl"}}
\description{
  This is the structure to hold the optimization control options.
  It takes the parameters, performs consistency checks, and returns the
  control structure.  Alternatively, it overwrites the supplied
  parameters in an existing \code{MaxControl} structure.  There is also
  a method to extract the control structure from the estimated
  \sQuote{maxim}-objects.  
}
\section{Slots}{
  The default values and definition of the slots:
  \describe{
    \item{\code{tol}}{1e-8, stopping condition 
      for \code{\link{maxNR}} and related optimizers.
      Stop if the absolute difference
      between successive iterations is less than \code{tol}, returns
      code 2.}
    \item{\code{reltol}}{sqrt(.Machine$double.eps), relative convergence
      tolerance (used by \code{\link{maxNR}} related optimizers, and
      \code{\link{optim}}-based optimizers.
      The algorithm stops if
      it is unable to increase the value by a factor of
      \code{reltol*(abs(val) + reltol)}
      at a step.  Returns code 2.}
    \item{\code{gradtol}}{1e-6, stopping condition
      for \code{\link{maxNR}} and related optimizers.
      Stop if norm of the gradient is
      less than \code{gradtol}, returns code 1.}
    \item{\code{steptol}}{1e-10, stopping/error condition
      for \code{\link{maxNR}} and related optimizers.
      If \code{qac == "stephalving"} and the quadratic
      approximation leads to a worse, instead of a better value, or to
      \code{NA}, the step length
      is halved and a new attempt is made.  If necessary, this procedure is repeated
      until step < \code{steptol}, thereafter code 3 is returned.}
    %
    \item{\code{lambdatol}}{1e-6, (for \code{\link{maxNR}} related
      optimizers)
      controls whether Hessian is treated as negative
      definite.  If the
      largest of the eigenvalues of the Hessian is larger than
      \code{-lambdatol} (Hessian is not negative definite),
      a suitable diagonal matrix is subtracted from the
      Hessian (quadratic hill-climbing) in order to enforce negative
      definiteness.}
    %
    \item{\code{qac}}{"stephalving", character, Qadratic Approximation
      Correction for \code{\link{maxNR}} related optimizers.  When the new
      guess is worse than the initial one, program attempts to correct it:
      \code{"stephalving"} decreases the
      step but keeps the direction.
      \code{"marquardt"} uses
      Marquardt (1963) method by decreasing the step length while also
      moving closer to the pure gradient direction.  It may be faster and
      more robust choice in areas where quadratic approximation behaves poorly.}
    \item{\code{qrtol}}{1e-10, QR-decomposition tolerance
      for Hessian inversion in \code{\link{maxNR}} related optimizers.
    }
    \item{\code{lambda0}}{1e-2, a positive numeric, initial correction term
      for \cite{Marquardt (1963)} correction in
      \code{\link{maxNR}}-related optimizers}
    \item{\code{lambdaStep}}{2, how much the \cite{Marquardt (1963)}
      correction is decreased/increased at
      successful/unsuccesful step
      for \code{\link{maxNR}} related optimizers}
    \item{\code{maxLambda}}{1e12, maximum allowed correction term
      for \code{\link{maxNR}} related optimizers.
      If exceeded, the
      algorithm exits with return code 3.}
    %
    \item{\code{alpha}}{1, Nelder-Mead simplex method reflection
      factor (see Nelder \& Mead, 1965)}
    \item{\code{beta}}{0.5, Nelder-Mead contraction factor}
    \item{\code{gamma}}{2, Nelder-Mead expansion factor}
    % SANN
    \item{\code{SANN_cand}}{\code{NULL} or a function for \code{"SANN"} algorithm
      to generate a new candidate point;
      if \code{NULL}, Gaussian Markov kernel is used
      (see argument \code{gr} of \code{\link{optim}}).}
    %
    \item{\code{iterlim}}{150, stopping condition.  Stop if more than \code{iterlim}
      iterations performed.  Note that \sQuote{iteration} may mean
      different things for different optimzers.}
    \item{\code{printLevel}}{0, the level of verbosity.  Larger values print
      more information.  Result depends on the optimizer.  Form
      \code{print.level} is also accepted by the methods for compatibility.}
  }
}
\section{Methods}{
  \describe{
    \item{maxControl}{\code{(\dots)} creates a \dQuote{MaxControl} object.  The
      arguments must be in the form \emph{option1 = value1, option2 =
	value2, ...}.}
    \item{maxControl}{\code{(x = "MaxControl", \dots)} overwrites parameters
      of and existing \dQuote{MaxControl} object}
    \item{maxControl}{\code{(x = "maxim")} extracts \dQuote{MaxControl}
      structure from an estimated model}
    \item{show}{shows the parameter values}
  }
}
\section{Details}{
  Typically, the control options are supplied in the form of a list, in which
  case the corresponding default values are overwritten by the
  user-specified ones.  However, one may also create the control
  structure by \code{maxControl(opt1=value1, opt2=value2, ...)} and
  supply such value directly to the optimizer.  In this case the
  optimization routine takes all the values from the control object.
  }
\references{
  \itemize{
    \item Nelder, J. A. & Mead, R. A (1965) Simplex Method for Function
    Minimization \emph{The Computer Journal} \bold{7}, 308--313
    \item Marquardt, D. W. (1963) An Algorithm for Least-Squares Estimation of
    Nonlinear Parameters \emph{Journal of the Society for Industrial and
    Applied Mathematics} \bold{11}, 431--441 
  }
}
\author{
  Ott Toomet \email{otoomet@gmail.com}
}
\note{
  Several control parameters can also be supplied directly to the
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
## create control object and use it for optimization
co <- maxControl(printLevel=2, qac="marquardt", lambda0=1)
res <- maxNR(quadForm, start=D, control=co)
print(summary(res))
## Now perform the same with no trace information
co <- maxControl(co, printLevel=0)
res <- maxNR(quadForm, start=D, control=co) # no tracing information
print(summary(res))  # should be the same as above
maxControl(res) # shows the control structure
}
\keyword{utilities}
