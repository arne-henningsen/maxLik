\name{sumt}
\Rdversion{1.1}
\alias{sumt}
\title{
  Equality-constrained optimization
}
\description{
  Sequentially Unconstrained Maximization Technique (SUMT) based
  optimization for linear equality constraints.

  This implementation is mostly intended to be called from other
  maximization routines, such as \code{\link{maxNR}}.
}
\usage{
sumt(fn, grad=NULL, hess=NULL,
start,
maxRoutine, constraints, SUMTTol = sqrt(.Machine$double.eps),
SUMTQ = 10,
SUMTRho0 = NULL,
print.level = 0, SUMTMaxIter = 100, ...)
}
\arguments{
  \item{fn}{
    function of a (single) vector parameter.  The function may have more
    arguments, but those are not treated as parameter
  }
  \item{grad}{
    gradient function of \code{fn}.  NULL if missing
  }
  \item{hess}{
    hessian matrix of the \code{fn}.  NULL if missing
  }
  \item{start}{
    initial value of the parameter.
}
  \item{maxRoutine}{
    maximization algorithm
}
  \item{constraints}{list, information for constrained maximization.
    Currently two components are supported: \code{eqA} and \code{eqB}
    for linear equality constraints: \eqn{A \beta + B = 0}{A \%*\% beta
      + B = 0}.  The user must ensure that the matrices \code{A} and
    \code{B} are conformable.}
  \item{SUMTTol}{
    stopping condition.  If the components of the parameter between
    successive iterations differ less
    than \code{SUMTTol}, the algorithm stops
  }
  \item{SUMTQ}{
    a double greater than one controlling the growth of the \code{rho}
    as described in Details. Defaults to 10.
  }
  \item{SUMTRho0}{
    Initial value for \code{rho}.  If not specified, a (possibly)
    suitable value is
    selected.  See Details.
    }
  \item{print.level}{
    Integer, debugging information.  Larger number print more details.
  }
  \item{SUMTMaxIter}{
    Maximum SUMT iterations
  }
  \item{\dots}{
    Other arguments to \code{maxRoutine} and \code{fn}.
  }
}
\details{
  The Sequential Unconstrained Minimization Technique is a heuristic
  for constrained optimization.  To minimize a function \eqn{f}{f}
  subject to
     constraints, one employs a non-negative function \eqn{P}{P} penalizing
     violations of the constraints, such that \eqn{P(x)}{P(x)} is zero
     iff \eqn{x}{x}
     satisfies the constraints.  One iteratively minimizes
     \eqn{L(x) + \varrho_k P(x)}{L(x) + rho_k P(x)}, where the
     \eqn{\varrho}{rho}
     values are increased according to the rule
     \eqn{\varrho_{k+1} = q \varrho_k}{rho_{k+1} = q rho_k} for some
     constant \eqn{q > 1}{q > 1}, until convergence is
     obtained in the sense that the Euclidean distance between
     successive solutions \eqn{x_k}{x_k} and \eqn{x_{k+1}}{x_{k+1}} is
     small enough.  Note that
     the "solution" \eqn{x}{x} obtained does not necessarily satisfy the
     constraints, i.e., has zero P(x).  Note also that there is no
     guarantee that global (approximately) constrained optima are
     found.  Standard practice would recommend to use the best solution
     found in "sufficiently many" replications of the algorithm.

     The unconstrained minimizations are carried out by either any of
     the maximization algorithms in the \pkg{maxLik}, such as
     \code{\link{maxNR}}.  Analytic gradient and hessian are used if
     provided, numeric ones otherwise.
}
\value{
  Object of class 'maxim'.  In addition, a component
  \item{constraints}{A list, describing the constrained optimization.
    Includes the following components:
    \itemize{
      \item{type}{type of constrained optimization}
      \item{outer.iterations}{number of iterations in the SUMT step}
      \item{barrier.value}{value of the penalty function at maximum}
    }
  }
}
\section{Note}{
  It may be a lot more efficient to embrace the actual function to be
  optimized to an outer function, which calculates the actual
  parameters based on a smaller set of parameters and the
  constraints.
}
\author{
  Ott Toomet \email{otoomet@ut.ee}
}
\seealso{
  \code{\link[pkg:clue]{sumt}}
}
\keyword{optimize}
