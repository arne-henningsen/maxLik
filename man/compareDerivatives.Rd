\name{compareDerivatives}
\alias{compareDerivatives}
\title{function to compare analytic and numeric derivatives}
\description{
  This function compares analytic and numerical derivative and prints related
  diagnostics information.  It is intended for testing and debugging
  code for analytic derivatives
  for maximisation algorithms.
}
\usage{
compareDerivatives(f, grad, hess=NULL, t0, eps=1e-6, print=TRUE, ...)
}
\arguments{
  \item{f}{
    function to be differentiated. The parameter (vector) of interest
    must be the first argument. The function may return a vector.
  }
  \item{grad}{
    analytic gradient.  This may be either a function,
    returning the analytic gradient, or a numeric vector, the pre-computed
    gradient.  The function must use the same set of
    parameters as \code{f}.  If \code{f} is a vector-valued function,
    grad must return/be a matrix where the number of rows equals the number
    of components of \code{f}, and the number of columns must equal to
    the number of components in \code{t0}.
  }
  \item{hess}{
    function returning the analytic hessian.  If present, hessian
    matrices are compared too.  Only appropriate for scalar-valued
    functions.
  }
  \item{t0}{
    parameter vector indicating the point at which the derivatives are
    compared. The derivative is taken with respect to this vector.
  }
  \item{eps}{
    numeric. Step size for numeric differentiation. Central derivative
    is used.
  }
  \item{print}{
    logical:  TRUE to print a summary, FALSE to return the comparison
    only (invisibly).
  }
  \item{\dots}{
    further arguments to \code{f}, \code{grad} and \code{hess}.
  }
}
\details{
  Analytic derivatives (and Hessian) substantially improve the
  estimation speed and reliability.  However, those are
  typically hard to program.  This utility compares the results and
  analytic and (internally calculated) numeric derivatives.  
  For every component of \code{f}, it prints the parameter value, analytic and
  numeric derivative, and their relative difference
  \deqn{\textrm{rel.diff} = \frac{\textrm{analytic} -
      \textrm{numeric}}{\frac{1}{2}(\textrm{analytic} + \textrm{numeric})}.}{rel.diff = (analytic -
    numeric)/(0.5*(analytic + numeric)).}
  If analytic = 0 = numeric, the rel.diff = 0.  If
  analytic derivatives are correct and the function is sufficiently
  smooth, expect the relative differences to be less than \eqn{10^{-7}}{1e-7}.
}
\value{
  A list with the following components:
  \item{t0}{the input argument \code{t0}}
  \item{f.t0}{f(t0)}
  \item{compareGrad}{
    a list with components \code{analytic} = grad(t0), \code{nmeric} =
    numericGradient(f, t0), and their \code{rel.diff}.
  }
  \item{maxRelDiffGrad}{max(abs(rel.diff))}

  If \code{hess} is also provided, the following optional components
  are also present:

  \item{compareHessian}{
    a list with components \code{analytic} = hess(t0), \code{numeric}
    = numericGradient(grad, t0), and their \code{rel.diff}.
  }
  \item{maxRelDiffHess}{max(abs(rel.diff)) for the Hessian}
}
\author{Ott Toomet \email{otoomet@ut.ee} and Spencer Graves}
\seealso{
  \code{\link{numericGradient}}
  \code{\link{deriv}}
}

\examples{
## A simple example with sin(x)' = cos(x)
f <- function(x)c(sin=sin(x))
Dsin <- compareDerivatives(f, cos, t0=c(angle=1))
D2sin <- compareDerivatives(f, cos, function(x)-sin(x), t0=1)

##
## Example of log-likelihood of normal density.  Two-parameter
## function.
##
x <- rnorm(100, 1, 2) # generate rnorm x
l <- function(b) sum(log(dnorm((x-b[1])/b[2])/b[2]))
              # b[1] = mu, b[2] = sigma
gradl <- function(b) {
    c(mu=sum(x - b[1])/b[2]^2,
    sigma=sum((x - b[1])^2/b[2]^3 - 1/b[2]))
}
gradl. <- compareDerivatives(l, gradl, t0=c(mu=1,sigma=2))

##
## An example with f returning a vector, t0 = a scalar
##
trig <- function(x)c(sin=sin(x), cos=cos(x))
Dtrig <- function(x)c(sin=cos(x), cos=-sin(x))
Dtrig. <- compareDerivatives(trig, Dtrig, t0=1)

D2trig <- function(x)-trig(x)
D2trig. <- compareDerivatives(trig, Dtrig, D2trig, t0=1)


}
\keyword{math}
\keyword{utilities}
