# _maxLik_: maximum likelihood estimation and related tools

This package is designed for Maximum Likelihood (ML) estimation.  It
contains a number of optimization routines, summary methods, and other
tools that are useful for ML estimation.

The highlights include:

* `maxLik`: the central function that performs the ML estimation.  It
  can be called with loglik function and start values as simply as
  `maxLik(loglik, start)`.
* support for _BHHH_ method
* tools to help debugging analytic gradient and Hessian
* a number of summary methods, including `summary`, `stdDev` and
  `tidy` for a quick summaries and standard deviations.


## History

The package (and its name) was inspired by the _maxlik_ library in
_gauss_ programming language.

The very first code of `maxLik` originates from a PhD econometrics
course in fall 2000.  The course was taught by Lars Muus at Aarhus
University, and a problem set asked the students to implement
Gauss-Newton method.  Later, OT could not understand error messages of
`nlm` function, and amended the Gauss-Newton to Newton-Raphson.  This
is the Newton-Raphson method that is one of the central optimizers in
current `maxLik`.
