maxBFGSYCCompute <- function(fn, grad=NULL, hess=NULL,
                         start, print.level=0,
                  tol=1e-6, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  iterlim=150,
                  fixed=NULL,
                  ...) {
   ##                       iterlim = 2000,
   ##                       tol = 1E-06,
   ## This function is originally developed by Yves Croissant (and placed in 'mlogit' package).
   ## Fitted for 'maxLik' by Ott Toomet
   ## 
   ## BFGS maximisation, implemented by Yves Croissant
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ##               If fn returns the value with attributes 'gradient'
   ##               and 'hessian', those are used instead of calculatin
   ##               new gradient and hessian values
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (numeric used if missing)
   ## start       - initial parameter vector (eventually w/names)
   ## steptol     - minimum step size
   ## ...         - extra arguments for fn()
   ## The stopping criteria
   ## tol         - maximum allowed absolute difference between sequential values
   ## reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ## gradtol     - maximum allowed norm of gradient vector
   ## iterlim     - maximum # of iterations
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ##               Other paramters are treated as variable (free).
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             activePar - logical vector, which parameters are active (not constant)
   ## activePar   logical vector, which parameters were treated as free (resp fixed)
   ## iterations  number of iterations
   ## type        "BFGS-YC maximisation"
   ## 
   func <- function(theta, sumObs = TRUE, ...) {
      f <- fn(theta, ...)
      if( sumObs ) {
         fAttr <- attributes(f)
         f <- sum(f)
         mostattributes(f) <- fAttr
      }
      return( f )
   }
   gradient <- function(theta, sumObs = TRUE,
                        suppliedValue=NULL, ...) {
      ## suppliedValue: use gradient value supplied from elsewhere
      ##           (attribute to fn) and only do sanity checks
      if(is.null(gr <- suppliedValue)) {
         if(!is.null(grad)) {  # use user-supplied if present
            gr <- grad(theta, ...)
         } else {
            gr <- numericGradient(f = func, t0 = theta,
                                  activePar=activePar, sumObs = sumObs, ...)
                           # Note we need nObs rows x nParam cols
         }
      }
      ## Now check if the gradient is vector or matrix...
      if( !is.null(dim(gr)) && sumObs ) {
         gr <- colSums(gr)
      } else {
         ## ... or vector if only one parameter
         if(length(gr) > nParam && sumObs ) {
            gr <- sum(gr)
         }
      }
      if( is.null( dim( gr ) ) ) {
         gr[ !activePar ] <- NA
      } else {
         gr[ , !activePar ] <- NA
      }
      return(gr)
   }
   ##
   maxim.type <- "BFGS-YC maximization"
   argNames <- c( "fn", "grad", "hess", "start", "print.level",
      "tol", "reltol", "gradtol", "steptol", 
      "iterlim", "activePar", "fixed" )
   checkFuncArgs( fn, argNames, "fn", "maxNR" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxNR" )
   }
  param <- start
   nimed <- names(start)
   nParam <- length(param)
   ## establish the active parameters.  Internally, we just use 'activePar'
   activePar <- !fixed
   rm( fixed )
   ##
  chi2 <- 1E+10
  iter <- 0
  # eval a first time the function, the gradient and the hessian
  x <- func(param, ...)
   if (print.level > 0)
    cat(paste("Initial value of the function :", as.numeric(x), "\n"))
   if(is.na(x)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   if(is.infinite(x) & (x > 0)) {
                                        # we stop at +Inf but not at -Inf
      result <- list(code=5, message=maximMessage("5"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   ##
   gri <- gradient(param, suppliedValue=attr(x, "gradient"), sumObs=FALSE, ...)
   gr <- colSums(gri)
                           # gradient by individual observations and simple
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(gr)
   }
   if(any(is.na(gr[activePar]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(gr[activePar]))) {
      stop("Infinite initial gradient")
   }
   if(length(gr) != nParam) {
      stop( "length of gradient (", length(gr),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   ##
   Hm1 <- solve(crossprod(gri[,activePar]))
                           # initial approximation of inverse Hessian (as in BHHH).
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(x), "\n")
      a <- cbind(start, gr, as.integer(activePar))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
   }
  direction <- rep(0, nParam)
   repeat {
      if( iter >= iterlim) {
         code <- 4; break
      }
    step <- 1
     direction[activePar] <- -as.vector(Hm1 %*% gr[activePar])
      iter <- iter + 1
    oldx <- x
     oldgr <- gr
    oldparam <- param
    param[activePar] <- oldparam[activePar] - step * direction[activePar]
    x <- func(param, ...)
    while((is.na(x) | x < oldx) & step > steptol) {
       step <- step/2
       if(print.level > 2) {
          cat("function values: old ", oldx, ", new ", x, ", difference ", x - oldx, " -> step ", step, "\n", sep="")
          if(print.level > 3) {
             resdet <- cbind(param = param, gradient = gr, direction=direction, active=activePar)
             print(resdet)
          }
       }
       param[activePar] <- oldparam[activePar] - step * direction[activePar]
       x <- func(param, ...)
    }
    if(step < steptol) {
                           # we did not find a better place to go...
       samm <- list(theta0=oldparam, f0=oldx, climb=direction)
    }
     gri <- gradient(param, suppliedValue=attr(x, "gradient"), sumObs=FALSE, ...)
     gr <- colSums(gri)
      incr <- step * direction
      y <- gr - oldgr
     Hm1 <- Hm1 +
        outer( incr[activePar], incr[activePar]) *
          (sum(y[activePar] * incr[activePar]) +
           as.vector( t(y[activePar]) %*% Hm1 %*% y[activePar])) / sum(incr[activePar] * y[activePar])^2 -
             (Hm1 %*% outer(y[activePar], incr[activePar])
              + outer(incr[activePar], y[activePar]) %*% Hm1)/
                  sum(incr[activePar] * y[activePar])
      chi2 <- -  crossprod(direction[activePar], oldgr[activePar])
    if (print.level > 0){
       cat("iteration ", iter, ", step = ",step,
                      ", lnL = ", x,", chi2 = ",
           chi2,"\n",sep="")
       if (print.level > 1){
          resdet <- cbind(param = param, gradient = gr, direction=direction, active=activePar)
          print(resdet)
          if(print.level > 3) {
             cat("Approximated Hessian:\n")
             print(Hm1)
          }
       }
       cat("--------------------------------------------\n")
    }
      if( step < steptol) {
         code <- 3; break
      }
      if( sqrt( t(gr[activePar])%*%gr[activePar]) < gradtol) {
         code <-1; break
      }
      if(x - oldx < tol) {
         code <- 2; break
      }
      if(x - oldx < reltol*(x + reltol)) {
         code <- 2; break
      }
      if(is.infinite(x) & x > 0) {
         code <- 5; break
      }
  }
   names(gr) <- colnames(gri) <- names(param)
  attr(x, "fixed") <- !activePar
  est.stat <- structure(list(elaps.time = NULL, nb.iter = iter, eps = chi2,
                             method = "BFGS-YC", message = message), class = 'est.stat')
  result <- list(optimum = x,
                 coefficients = param,
                 est.stat = est.stat
                 )
  result
}


maxBFGSYC <- function(fn, grad=NULL, hess=NULL, start, print.level=0,
                  tol=1e-8, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  iterlim=150,
                  constraints=NULL,
                  fixed=NULL,
                  activePar=NULL,
                  ...) {
   ## Newton-Raphson maximisation
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (numeric used if missing)
   ## start       - initial parameter vector (eventually w/names)
   ## steptol     - minimum step size
   ## ...         - extra arguments for fn()
   ## The stopping criteria
   ## tol         - maximum allowed absolute difference between sequential values
   ## reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ## gradtol     - maximum allowed norm of gradient vector
   ## iterlim     - maximum # of iterations
   ## activePar   - an index vector -- which parameters are taken as
   ##               variable (free).  Other paramters are treated as
   ##               fixed constants
   ## fixed         index vector, which parameters to keep fixed
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             activePar - logical vector, which parameters are active (not constant)
   ## activePar   logical vector, which parameters were treated as free (resp fixed)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   argNames <- c( "fn", "grad", "hess", "start", "print.level",
      "tol", "reltol", "gradtol", "steptol",
      "iterlim", "activePar", "fixed" )
   checkFuncArgs( fn, argNames, "fn", "maxNR" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxNR" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxNR" )
   }
   ## establish the active parameters.  Internally, we just use 'activePar'
   fixed <- prepareFixed( start = start, activePar = activePar,
      fixed = fixed )
   if(is.null(constraints)) {
       result <- maxBFGSYCCompute(fn=fn, grad=grad, hess=hess,
                              start=start,
                              print.level=print.level,
                              tol=tol, reltol=reltol,
                              gradtol=gradtol, steptol=steptol,
                              iterlim=iterlim,
                              fixed=fixed,
                              ...)
   } else {
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         stop("Inequality constraints not implemented for maxNR")
      } else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=fn, grad=grad, hess=hess,
                        start=start,
                        maxRoutine=maxBFGSYC,
                        constraints=constraints,
                        print.level=print.level,
                        tol=tol, reltol=reltol,
                        gradtol=gradtol, steptol=steptol,
                        iterlim=iterlim,
                        fixed=fixed,
                        ...) 
      } else {
         stop("maxBFGS only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }
   return( result )
}
