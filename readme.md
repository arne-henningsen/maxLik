# _maxLik_: maximum likelihood estimation and related tools

This package contains a set of functions and tools for Maximum
Likelihood (ML) estimation.  The focus of the package is on non-linear
optimization from the ML viewpoint, and it provides several
convenience wrappers and tools, like _BHHH_ algorithm,
variance-covariance matrix, standard errors, and summary methods.

## Highlights:

* `maxLik`: the central function that performs the ML estimation.  It
  can be called with loglik function and start values as simply as
  `maxLik(loglik, start)`.
* a number of optimization methods with unified interface, most of
  which can be called from through `maxLik`.
* support for _BHHH_ method
* tools to help debugging analytic gradient and Hessian
* a number of summary methods, including `summary`, `stdDev` and
  `tidy` for a quick summaries and standard deviations.

## Details

_maxLik_ package is a set of convenience tools and wrappers
focusing on
Maximum Likelihood (ML) analysis, but it also contains tools for
other optimization tasks.
The package includes a) wrappers for several
existing optimizers (implemented by _stats::optim_; b) original
optimizers, including Newton-Raphson and Stochastic Gradient Ascent;
and c) several convenience tools
to use these optimizers from the ML perspective.  Examples are _BHHH_
optimization (_maxBHHH_) and utilities that extract
standard errors from the estimates.  Other highlights include a unified
interface for all included optimizers, tools to test user-provided analytic
derivatives, and constrained optimization.

A good starting point to learn about the usage of _maxLik_ are the
included vignettes "Introduction: what is maximum likelihood",
"Maximum likelihood estimation with maxLik" and
"Stochastic Gradient Ascent in maxLik".  Another good
source is
Henningsen & Toomet (2011), an introductory paper to the package.
Use
`vignette(package="maxLik")` to see the available vignettes, and
`vignette("using-maxlik")` to read the usage vignette.

## Quick intro

From the user's perspective, the
central function in the package is `maxLik`.  In its
simplest form it takes two arguments: the log-likelihood function, and
a vector of initial parameter values (see the example below).
It returns an object of class
_maxLik_ with convenient methods such as
_summary_, _coef_, _stdEr_.  It also supports a plethora
of other arguments, for instance one can supply analytic gradient and
Hessian, select the desired optimizer, and control the optimization in
different ways.

A useful utility functions in the package is
`compareDerivatives` that
allows one to compare the analytic and numeric derivatives for debugging
purposes.
Another useful function is `condiNumber` for
analyzing multicollinearity problems in the estimated models.

In the interest of providing a unified user interface, all the
optimizers are implemented as maximizers in this package.  This includes
the _optim_-based methods, such as _maxBFGS_ and
_maxSGA_, the maximizer version of popular Stochastic
Gradient Descent.

An quick example:
```r
### estimate mean and variance of normal random vector

## create random numbers where mu=0, sd=1
set.seed(123)
x <- rnorm(50)

### log likelihood function.
llf <- function(param) {
   ## Note: 'param' is a 2-vector c(mu, sd)
   mu <- param[1]
   sd <- param[2]
   llValue <- dnorm(x, mean=mu, sd=sd, log=TRUE)
   sum(llValue)
}

## Estimate it with mu=1, sd=2 as start values
ml <- maxLik(llf, start = c(mu=1, sigma=2) )
print(summary(ml))
## Estimates close to c(0, 1) :-)
```

## Authors

* Ott Toomet
* Arne Henningsen
* Spencer Graves
* Yves Croissant
* David Hugh-Jones
* Lucca Scrucca

Maintainer: Ott Toomet

## References

Henningsen A, Toomet O (2011). _maxLik: A package for maximum
    likelihood estimation in R._
  Computational Statistics, 26(3), 443-458. doi:
10.1007/s00180-010-0217-1.



## History

The package (and its name) was inspired by the _maxlik_ library in
_gauss_ programming language.

The very first code of `maxLik` originates from a PhD econometrics
course in fall 2000.  The course was taught by Lars Muus at Aarhus
University, and a problem set asked the students to implement
Gauss-Newton method.  Later, when OT could not understand an error message of
`nlm` function, he converted the Gauss-Newton function to Newton-Raphson.  This
is the Newton-Raphson method that is one of the central optimizers in
current `maxLik`.
