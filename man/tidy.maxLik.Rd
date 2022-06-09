\name{tidy.maxLik}
\alias{tidy.maxim}
\alias{tidy.maxLik}
\alias{glance.maxim}
\alias{glance.maxLik}
\title{tidy and glance methods for maxLik objects}
\description{
  These methods return quick
  summary information about the estimated model, either one parameter
  per row (tibble), or on a single row (glance).
  Both require the \pkg{tibble} package to be installed.
}
\usage{
\method{tidy}{maxim}(x,  ...)
\method{tidy}{maxLik}(x,  ...)
\method{glance}{maxim}(x, ...)
\method{glance}{maxLik}(x, ...)
}
\arguments{
  \item{x}{
    object of class 'maxLik' or 'maxim'.
  }
  \item{\ldots}{Not used.}
}
\value{
  \code{tidy()} for 'maxLik' objects returns a tibble
  (data frame) with one row per
  parameter, and with the following columns:
  \describe{
    \item{term}{The name of the estimated parameter (parameters are
      sequentially numbered if names missing).}
    \item{estimate}{The estimated parameter.}
    \item{std.error}{The standard error of the estimate.}
    \item{statistic}{The \eqn{z}{z}-statistic of the estimate.}
    \item{p.value}{The \eqn{p}{p}-value.}
  }
  This is essentially the same table as \code{summary}-method prints,
  just in form of a tibble (data frame).

  For 'maxim' objects, the corresponding tibble contains
  \describe{
    \item{term}{The name of the estimated parameter (parameters are
      sequentially numbered if names missing).}
    \item{estimate}{The estimated parameter.}
    \item{gradient}{The corresponding gradient value.}
  }
  
  \code{glance()} for 'maxLik' objects returns an
  one-row tibble with columns:
    \describe{
      \item{df}{The degrees of freedom of the model.}
      \item{logLik}{The log-likelihood of the model.}
      \item{AIC}{Akaike's Information Criterion for the model.}
      \item{nobs}{The number of observations, if this is available, otherwise \code{NA}.}
    }

    For 'maxim' objects the tibble contains
    \describe{
      \item{value}{The function value at optimum.}
      \item{iterations}{The number of iterations spent (what exactly is
	\dQuote{iteration} depends on the optimizer).
      }
      \item{gradient}{
	Norm of gradient at the optimum.
      }
      \item{code}{
	optimizer return code (see \code{\link{returnCode}}).
      }
    }
}
\seealso{
  The functions \code{\link[generics:tidy]{tidy}} and
  \code{\link[generics:glance]{glance}} in package \pkg{generics}, and
    \code{\link[=summary.maxLik]{summary}} to display the
    \dQuote{standard} summary information.
}
\author{David Hugh-Jones}
\examples{
## Example with a single parameter
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
a <- maxLik(loglik, start=2)
tidy(a)
glance(a)
## Example with a parameter vector
x <- rnorm(100)
loglik <- function(theta) {
   dnorm(x, mean=theta[1], sd=theta[2], log=TRUE)
}
a <- maxLik(loglik, start=c(mu=0, sd=1))
tidy(a)
glance(a)
}
