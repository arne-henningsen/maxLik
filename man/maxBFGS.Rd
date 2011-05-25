\name{maxBFGS}
\alias{maxBFGS}
\alias{maxSANN}
\alias{maxNM}

\title{BFGS, SANN and Nelder-Mead Maximization}
\description{
  These functions are wrappers for \code{\link{optim}} where the arguments are
  compatible with \code{\link{maxNR}}.  Note that there is a
  \code{\link{maxNR}}-based BFGS implementation \code{\link{maxBFGSR}}.
}
\usage{
maxBFGS(fn, grad = NULL, hess=NULL, start, fixed = NULL,
   print.level = 0, iterlim = 200, constraints = NULL,
   tol = 1e-08, reltol=tol,
   finalHessian=TRUE,
   parscale=rep(1, length=length(start)), ... )

maxSANN(fn, grad = NULL, hess = NULL, start, fixed = NULL,
   print.level = 0, iterlim = 10000, constraints = NULL,
   tol = 1e-08, reltol=tol,
   finalHessian=TRUE,
   cand = NULL, temp = 10, tmax = 10,
   parscale = rep(1, length = length(start)),
   random.seed = 123, ... )

maxNM(fn, grad = NULL, hess = NULL, start, fixed = NULL,
   print.level = 0, iterlim = 500, constraints = NULL,
   tol = 1e-08, reltol=tol,
   finalHessian=TRUE,
   parscale = rep(1, length = length(start)),
   alpha = 1, beta = 0.5, gamma = 2, ...)
}

\arguments{
  \item{fn}{function to be maximised.  Must have the parameter vector as
    the first argument.  In order to use numeric gradient
    and BHHH method, \code{fn} must return vector of
    observation-specific likelihood values.  Those are summed by maxNR
    if necessary.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of the function.  Must have the parameter vector as
    the first argument.  If \code{NULL}, numeric
    gradient is used (only maxBFGS uses gradient).  Gradient may return
    a matrix, where columns correspond to the parameters and rows to the
    observations (useful for maxBHHH).  The columns are summed internally.}
  \item{hess}{Hessian of the function.  Not used by any of these methods, for
    compatibility with \code{\link{maxNR}}.}
  \item{start}{initial values for the parameters.}
  \item{fixed}{parameters that should be fixed at their starting values:
      either a logical vector of the same length as argument \code{start},
      a numeric (index) vector indicating the positions of the fixed parameters,
      or a vector of character strings indicating the names of the fixed
      parameters (parameter names are taken from argument \code{start}).}
  \item{print.level}{a larger number prints more working information.}
  \item{iterlim}{maximum number of iterations.}
  \item{constraints}{either \code{NULL} for unconstrained optimization
    or a list with two components.  The components may be either
    \code{eqA} and \code{eqB} for equality-constrained optimization
    \eqn{A \theta + B = 0}{A \%*\% theta + B = 0}; or \code{ineqA} and
    \code{ineqB} for inequality constraints \eqn{A \theta + B \ge 0}{A
      \%*\% theta + B >= 0}.  The equality-constrained problem is forwarded
    to \code{\link{sumt}}, the inequality-constrained case to
    \code{\link{constrOptim2}}.
  }
  \item{tol, reltol}{the relative convergence tolerance (see
    \code{\link{optim}}).  \code{tol} is for compatibility with
    \code{\link{maxNR}}.}
  \item{finalHessian}{how (and if) to calculate the final Hessian.  Either
  \code{FALSE} (not calculate), \code{TRUE} (use analytic/numeric
  Hessian) or \code{"bhhh"}/\code{"BHHH"} for information equality
  approach.  The latter approach is only suitable for maximizing
  log-likelihood function.  It requires the gradient/log-likelihood to
  be supplied by individual observations, see \code{\link{maxBHHH}} for
  details. 
 }
  \item{cand}{a function used in the \code{"SANN"} algorithm
    to generate a new candidate point;
    if it is \code{NULL}, a default Gaussian Markov kernel is used
    (see argument \code{gr} of \code{\link{optim}}).}
  \item{temp}{controls the '"SANN"' method. It is the starting
          temperature for the cooling schedule. Defaults to '10'.}
  \item{tmax}{is the number of function evaluations at each temperature
          for the '"SANN"' method. Defaults to '10'. (see
    \code{\link{optim}})}
  \item{random.seed}{an integer used to seed R's random number generator.
    This is to ensure replicability when the \sQuote{Simulated Annealing}
    method is used. Defaults to 123.}
  \item{parscale}{A vector of scaling values for the parameters.
          Optimization is performed on 'par/parscale' and these should
          be comparable in the sense that a unit change in any element
          produces about a unit change in the scaled value. (see
    \code{\link{optim}})}
  \item{alpha, beta, gamma}{Scaling parameters for the
          '"Nelder-Mead"' method. 'alpha' is the reflection factor
          (default 1.0), 'beta' the contraction factor (0.5) and
          'gamma' the expansion factor (2.0). (see
    \code{\link{optim}})}
  \item{\dots}{further arguments for \code{fn} and \code{grad}.}
}

\details{
   The \sQuote{state} (or \sQuote{seed}) of R's random number generator
   is saved at the beginning of the \code{maxSANN} function
   and restored at the end of this function
   so that this function does \emph{not} affect the generation of random numbers
   although the random seed is set to argument \code{random.seed}
   and the \sQuote{SANN} algorithm uses random numbers.
}

\value{
  Object of class "maxim":
  \item{maximum}{value of \code{fn} at maximum.}
  \item{estimate}{best set of parameters found.}
  \item{gradient}{vector, gradient at parameter value \code{estimate}.}
  \item{gradientObs}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{value of Hessian at optimum.}
  \item{code}{integer. Success code, 0 is success (see
    \code{\link{optim}}).}
  \item{message}{character string giving any additional information
    returned by the optimizer, or NULL.}
  \item{fixed}{logical vector indicating which parameters are treated
    as constants.}
  \item{iterations}{two-element integer vector giving the number of
    calls to \code{fn} and \code{gr}, respectively.
    This excludes those calls needed to
    compute the Hessian, if requested, and any calls to \code{fn} to compute a
    finite-difference approximation to the gradient.}
  \item{type}{character string "BFGS maximisation".}
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \itemize{
      \item{type}{type of constrained optimization}
      \item{outer.iterations}{number of iterations in the constraints step}
      \item{barrier.value}{value of the barrier function}
    }
  }
}
\author{Ott Toomet \email{otoomet@ut.ee}, Arne Henningsen}

\seealso{\code{\link{optim}}, \code{\link{nlm}}, \code{\link{maxNR}},
  \code{\link{maxBHHH}}, \code{\link{maxBFGSR}}.}

\examples{
# Maximum Likelihood estimation of the parameter of Poissonian distribution
n <- rpois(100, 3)
loglik <- function(l) n*log(l) - l - lfactorial(n)
# we use numeric gradient
summary(maxBFGS(loglik, start=1))
# you would probably prefer mean(n) instead of that ;-)
# Note also that maxLik is better suited for Maximum Likelihood
###
### Now an example of constrained optimization
###
f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may want to use exp(- theta \%*\% theta) instead ;-)
}
## use constraints: x + y >= 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- maxNM(f, start=c(1,1), constraints=list(ineqA=A, ineqB=B),
print.level=1)
print(summary(res))
}
\keyword{optimize}
