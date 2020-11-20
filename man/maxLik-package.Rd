\name{maxLik-package}
\alias{maxLik-package}
\docType{package}
\title{
Maximum Likelihood Estimation
}
\description{
  This package contains a set of functions and tools for Maximum Likelihood (ML)
  estimation.  The focus of the package is on the non-linear
  optimization from the ML viewpoint, and it provides several convenience wrappers
  and tools, like BHHH algorithm and extraction of variance-covariance
  matrix.
}
\details{
\pkg{maxLik} package is a set of convenience tools and wrappers
focusing on
Maximum Likelihood (ML) analysis, but it also contains tools for
other optimization tasks.
The package includes a) wrappers for several
existing optimizers (implemented by \code{\link{optim}}); b) original
optimizers, including Newton-Raphson and Stochastic Gradient Ascent;
and c) several convenience tools
to use these optimizers from the ML perspective.  Examples are BHHH
optimization (\code{\link{maxBHHH}}) and utilities that extract
standard errors from the estimates.  Other highlights include a unified
interface for all included optimizers, tools to check the programmed analytic
derivatives, and constrained optimization.

A good starting point to learn about the usage of \pkg{maxLik} are
Henningsen \& Toomet (2011), and the provided vignette
\dQuote{Stochastic Gradient Ascent in maxLik}.

From the user's perspective, the
central function in the package is \code{\link{maxLik}}.  In the
simplest form it takes two arguments: the log-likelihood function, and
a vector of initial parameter values.  It returns an object of class
\sQuote{maxLik} with convenient methods such as
\code{\link[maxLik:summary.maxLik]{summary}},
\code{\link[maxLik:maxLik-methods]{coef}}, and
\code{\link[maxLik:maxLik-methods]{stdEr}}.  It also supports a plethora
of other arguments, for instance one can supply analytic gradient and
Hessian, select the desired optimizer, and control the optimization in
different ways.

One of the most useful utility functions in the package is
\code{\link{compareDerivatives}} that
allows one to compare the analytic and numeric derivatives for debugging
purposes.
Another useful function is \code{\link{condiNumber}} for
analyzing multicollinearity problems in the estimated models.

In the interest of providing a unified user interface, all the
optimizers are implemented as maximizers in this package.  This includes
the \code{\link{optim}}-based methods, such as \code{\link{maxBFGS}} and
\code{\link{maxSGA}}, the maximizer version of popular Stochastic
Gradient Descent.
}
\author{
Ott Toomet <otoomet@gmail.com>,
Arne Henningsen <arne.henningsen@gmail.com>,
with contributions from Spencer Graves and Yves Croissant

Maintainer: Ott Toomet <otoomet@gmail.com>
}
\references{
  Henningsen A, Toomet O (2011). \dQuote{maxLik: A package for maximum
    likelihood estimation in R.}
  Computational Statistics, 26(3), 443-458. doi:
  \href{http://dx.doi.org/10.1007/s00180-010-0217-1}{10.1007/s00180-010-0217-1}.
}
\keyword{Basics|package}
\keyword{Mathematics|optimize}
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
## estimate mean and variance of normal random vector
set.seed(123)
x <- rnorm(50, 1, 2 )

## log likelihood function.
## Note: 'param' is a vector
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   llValue <- dnorm(x, mean=mu, sd=sigma, log=TRUE)
   return(sum(llValue))
}

## Estimate it.  Take standard normal as start values
ml <- maxLik( llf, start = c(mu=0, sigma=1) )
print(summary(ml))
## Estimates close to c(1,2) :-)

## Example how to use maxLik in your own function and allow users
## to override the default parameters
##
## 'estimate': user contructed estimation routine
## Note: it accepts both 'control' and '...'
estimate <- function(control=NULL, ...) {
   res <- maxLik(llf, start=c(mu=1, sigma=1),
                 control=c(list(iterlim=100), control),
                           # user-supplied 'control' overrides default
                           # 'iterlim=100'
                 ...)
   res			  
}
m <- estimate(control=list(iterlim=1), fixed=2)
                           # user can override default 'iterlim' and
                           # supply additional parameters ('fixed')
show(maxControl(m))
                           # iterlim should be 1
print(coef(m))
                           # sigma should be 1.000

}
