\name{maxNR}
\alias{maxNR}
\title{Newton-Raphson Maximization}
\description{
  Unconstrained maximization based on Newton-Raphson method.
}
\usage{
maxNR(fn, grad = NULL, hess = NULL, start, print.level = 0,
      tol = 1e-08, reltol=sqrt(.Machine$double.eps), gradtol = 1e-06,
      steptol = 1e-10, lambdatol = 1e-06, qrtol = 1e-10,
      iterlim = 150,
      constraints=NULL,
      fixed=NULL, activePar=NULL, ... )
}
\arguments{
  \item{fn}{function to be maximized.  It must return either a single
    number or a numeric vector, which is summed (this is useful for BHHH
    method).  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.

    \code{fn} may also return attributes "gradient" and/or "hessian".
    If these attributes are set, the algorithm uses the corresponding
    values as
    gradient and hessian.  
  }
  \item{grad}{gradient of the function.  If \code{NULL}, numeric
    gradient is used.  It must return a gradient vector, or matrix where
    \emph{columns} correspond to individual parameters.  The column sums
    are treated as gradient components (this is useful for BHHH method).
  }
  \item{hess}{Hessian matrix of the function.  If missing, numeric
    Hessian, based on gradient, is used.  }
  \item{start}{initial value for the parameter vector}
  \item{print.level}{this argument determines the level of
    printing which is done during the minimization process. The default
    value of 0 means that no printing occurs, a value of 1 means that
    initial and final details are printed and a value of 2 means that
    full tracing information for every iteration is printed.  Higher
    values will result in even more details.
  }
  \item{tol}{stopping condition.  Stop if the absolute difference
    between successive iterations less than \code{tol}, return
    \code{code=2}.}
  \item{reltol}{Relative convergence tolerance.  The algorithm stops if
          it is unable to increase the value by a factor of 'reltol *
          (abs(val) + reltol)' at a step.  Defaults to
          'sqrt(.Machine\$double.eps)', typically about '1e-8'.}
  \item{gradtol}{stopping condition.  Stop if the norm of the gradient
    less than \code{gradtol}, return \code{code=1}.}
  \item{steptol}{stopping/error condition.  If the quadratic
    approximation leads to lower function value instead of higher, or
    \code{NA}, the step length
    is halved and a new attempt is made.  This procedure is repeated
    until step < \code{steptol}, thereafter \code{code=3} is returned.}
  \item{lambdatol}{control whether the Hessian is treated as negative
    definite.  If the
    largest of the eigenvalues of the Hessian is larger than
    \code{-lambdatol}, a suitable diagonal matrix is subtracted from the
    Hessian (quadratic hill-climbing) in order to enforce nagetive
    definiteness.}
  \item{qrtol}{QR-decomposition tolerance}
  \item{iterlim}{stopping condition.  Stop if more than \code{iterlim}
    iterations, return \code{code=4}.}
  \item{constraints}{either \code{NULL} for unconstrained optimization
    or a list with two components
    \code{eqA} and \code{eqB} for equality-constrained optimization
    \eqn{A \theta + B = 0}{A \%*\% theta + B = 0}.  The constrained
    problem is forwarded
    to \code{\link{sumt}}.
  }
  \item{fixed}{parameters that should be fixed at their starting values:
      either a logical vector of the same length as argument \code{start},
      a numeric (index) vector indicating the positions of the fixed parameters,
      or a vector of character strings indicating the names of the fixed
      parameters (parameter names are taken from argument \code{start}).}
  \item{activePar}{this argument is retained for backward compatibility only;
      please use argument \code{fixed} instead.}
  \item{\dots}{further arguments to \code{fn}, \code{grad} and
    \code{hess}.}
}
\details{
The algorithm can work with constant parameters and related changes of
parameter values.  Constant parameters are useful if a parameter value
is converging toward the boundary of support, or for testing.  

One way is to put
\code{activePar} to non-NULL, specifying which parameters we want to
change.  However, parameters can also be fixed in runtime by by
signaling with
\code{fn}.  This may be useful if an estimation converges toward the
edge of the parameter space possibly causing problems.  The value of
\code{fn} may have following attributes:
\itemize{
\item{constPar}{index vector.  Which parameters are redefined to
  constant}
\item{newVal}{a list with following components:
  \itemize{
  \item{index}{which parameters will have a new value}
  \item{val}{the new value of parameters}
  }
}
}
Hence, \code{constVal} specifies which parameters are treated as
constants.  \code{newVal} allows one to overwrite the existing parameter
values, possibly the non-constant values as well.  If the attribute
\code{newVal} is present, the new function value need not to exceed 
the previous one (maximization is not performed in that step).
}
\value{
  list of class "maxim" with following components:
  \item{maximum}{\code{fn} value at maximum (the last calculated value
    if not converged).}
  \item{estimate}{estimated parameter value.}
  \item{gradient}{last gradient value which was calculated.  Should be
    close to 0 if normal convergence.}
  \item{gradientObs}{matrix of gradients at parameter value \code{estimate}
      evaluated at each observation (only if \code{grad} returns a matrix
      or \code{grad} is not specified and \code{fn} returns a vector).}
  \item{hessian}{Hessian at the maximum (the last calculated value if
    not converged).}
  \item{code}{return code:
    \itemize{
    \item{1}{gradient close to zero (normal convergence).}
    \item{2}{successive function values within tolerance limit (normal
      convergence).}
    \item{3}{last step could not find higher value (probably not
      converged).}
    \item{4}{iteration limit exceeded.}
    \item{5}{Infinite value.}
    \item{6}{Infinite gradient.}
    \item{7}{Infinite Hessian.}
    \item{100}{Initial value out of range.}
    }
  }
  \item{message}{a short message, describing \code{code}.}
  \item{last.step}{list describing the last unsuccessful step if
    \code{code=3} with following components:
    \itemize{
    \item{theta0}{previous parameter value}
    \item{f0}{\code{fn} value at \code{theta0}}
    \item{climb}{the movement vector to the maximum of the quadratic approximation}
    }
  }
  \item{activePar}{logical vector, which parameters are not constants.}
  \item{iterations}{number of iterations.}
  \item{type}{character string, type of maximization.}
  \item{constraints}{A list, describing the constrained optimization
    (\code{NULL} if unconstrained).  Includes the following components:
    \itemize{
      \item{type}{type of constrained optimization}
      \item{outer.iterations}{number of iterations in the constraints step}
      \item{barrier.value}{value of the barrier function}
    }
  }
}
\references{Greene,  W., 2002, \emph{Advanced Econometrics}
  Goldfeld, S. M. and Quandt, R. E., 1972,
  \emph{Nonlinear Methods in Econometrics}.  Amsterdam: North-Holland}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{nlm}} for Newton-Raphson optimization,
  \code{\link{optim}} for different gradient-based optimization
  methods.}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta) sum(1/theta - t)
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and Hessian
a <- maxNR(loglik, start=1, print.level=2)
summary(a)
## You would probably prefer 1/mean(t) instead ;-)
## Estimate with analytic gradient and Hessian
a <- maxNR(loglik, gradlik, hesslik, start=1)
summary(a)
##
## Next, we give an example with vector argument:  Estimate the mean and
## variance of a random normal sample by maximum likelihood
## Note: you might want to use maxLik instead
##
loglik <- function(param) {
  mu <- param[1]
  sigma <- param[2]
  ll <- -0.5*N*log(2*pi) - N*log(sigma) - sum(0.5*(x - mu)^2/sigma^2)
  ll
}
x <- rnorm(1000, 1, 2) # use mean=1, stdd=2
N <- length(x)
res <- maxNR(loglik, start=c(0,1)) # use 'wrong' start values
summary(res)
###
### Now an example of constrained optimization
###
## We maximize exp(-x^2 - y^2) where x+y = 1
f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  exp(-(x^2 + y^2))
  ## Note: you may want to use exp(- theta \%*\% theta) instead ;-)
}
## use constraints: x + y = 1
A <- matrix(c(1, 1), 1, 2)
B <- -1
res <- maxNR(f, start=c(0,0), constraints=list(eqA=A, eqB=B), print.level=1)
print(summary(res))
}
\keyword{optimize}
