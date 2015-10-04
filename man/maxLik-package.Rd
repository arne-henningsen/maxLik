\name{maxLik-package}
\alias{maxLik-package}
\docType{package}
\title{
Maximum Likelihood Estimation
}
\description{
  This is a set of functions and tools to perform Maximum Likelihood (ML)
  estimation.  The focus of the package is on the non-linear
  optimization from the ML viewpoint, and it provides several convenience wrappers
  and tools, like BHHH algorithm and extraction of variance-covariance
  matrix.
}
\details{
\dQuote{maxLik} package is a set of convenience tools and wrappers to perform
Maximum Likelihood (ML) analysis.  It includes a) wrappers for several
existing optimizers (implemented by \code{\link{optim}}); b) original
optimizers, including Newton-Raphson; and c) several convenience tools
to use these optimizers from the ML perspective.  Examples are BHHH
optimization (\code{\link{maxBHHH}}) and utilities that extract
standard errors from the estimates.  Other highlights include a unified
interface for all included optimizers, tools to check the programmed analytic
derivatives, and constrained optimization.

From the user's perspective, the
central function in the package is \code{\link{maxLik}}.  In the
simplest form it takes two arguments: the log-likelihood function, and
a vector of parameters' start values.  It returns an object of class
\sQuote{maxLik} with convenient methods such as
\code{\link[maxLik:summary.maxLik]{summary}},
\code{\link[maxLik:coef.maxim]{coef}}, and
\code{\link[maxLik:stdEr.maxLik]{stdEr}}.  It also supports a plethora
of other arguments, for instance one can supply analytic gradient and
Hessian, select the desired optimizer, and control the optimization in
different ways.

One of the most useful utility functions in the package is
\code{\link{compareDerivatives}} that
allows one to compare the analytic and numeric derivatives for debugging
the derivative code.  Another useful function is \code{\link{condiNumber}} for
analyzing multicollinearity problems in the estimated models.
}
\author{
Ott Toomet <otoomet@gmail.com>,
Arne Henningsen <arne.henningsen@gmail.com>,
with contributions from Spencer Graves and Yves Croissant

Maintainer: Ott Toomet <otoomet@gmail.com>
}
%% \references{
%% }
\keyword{Basics|package }
\keyword{Mathematics|optimize}
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
## estimate mean and variance of normal random vector
set.seed( 123 )
x <- rnorm(50, 1, 2 )

## log likelihood function.
## Note: 'param' is a vector
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   N <- length( x )
   llValue <- -0.5 * N * log( 2 * pi ) - N * log( sigma ) -
      0.5 * sum( ( x - mu )^2 / sigma^2 )
   return(sum(llValue))
}

## Estimate it.  Take standard normal as start values
ml <- maxLik( llf, start = c(mu=0, sigma=1) )
print(summary(ml))
## Estimates close to c(1,2) :-)
}
