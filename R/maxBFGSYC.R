mlogit.optimX <- function(f, start,
                         method = c('bfgs', 'nr', 'bhhh'),
                         iterlim = 2000,
                         tol = 1E-06,
                         print.level = 1,
                         constPar = NULL,
                         ...){
  # construct a call for the function
  param <- start
  method <- match.arg(method)
  callT <- match.call(expand.dots = TRUE)
  f <- callT
  optimoptions <- c('iterlim', 'tol', 'method', 'print.level', 'constPar')
  m <- match(optimoptions, names(callT), 0L)
  if (sum(m)) f <- f[-m]
  f[[1]] <- as.name(f[[2]])
  K <- length(param)
  fixed <- rep(FALSE, K)
  if (!is.null(constPar)) fixed[constPar] <- TRUE
  f$gradient <- TRUE
  if (method == 'nr') f$hessian <- TRUE else f$hessian <- FALSE
  f[[2]] <- NULL
  names(f)[2] <- 'param'
  chi2 <- 1E+10
  i <- 0
  # eval a first time the function, the gradient and the hessian
  x <- eval(f, parent.frame())
  if (print.level > 0)
    cat(paste("Initial value of the function :", as.numeric(x), "\n"))
  g <- attr(x, "gradient")
  Hm1 <- solve(crossprod(attr(x, "gradi")[, !fixed]))
  d <- rep(0, K)
  while(abs(chi2) > tol && i < iterlim){
    if (method == "bfgs") d[!fixed] <- - as.vector(Hm1 %*% g[!fixed])
    else d[!fixed] <- - as.vector(solve(H, g[!fixed]))
    i <- i + 1
    oldx <- x
    oldg <- g
    f$param <- param
    f$direction <- d
    f$initial.value <- x
    x <- eval(f, parent.frame())
    g <- attr(x, "gradient")
    step <- attr(x, "step")
    param[!fixed] <- param[!fixed] + step * d[!fixed]
    if (method == 'bfgs'){
      incr <- step * d
      y <- g - oldg
      Hm1 <- Hm1 +
        outer( incr[!fixed], incr[!fixed]) *
          (sum(y[!fixed] * incr[!fixed]) +
           as.vector( t(y[!fixed]) %*% Hm1 %*% y[!fixed])) / sum(incr[!fixed] * y[!fixed])^2 -
             (Hm1 %*% outer(y[!fixed], incr[!fixed])
              + outer(incr[!fixed], y[!fixed]) %*% Hm1)/
                sum(incr[!fixed] * y[!fixed])
    }
    chi2 <- -  crossprod(d[!fixed], oldg[!fixed])
    if (print.level > 0){
      chaine <- paste("iteration ",i,", step = ",step,
                      ", lnL = ",round(x,8),", chi2 = ",
                      round(chi2,8),"\n",sep="")
    cat(chaine)
    }
    if (print.level > 1){
      resdet <- rbind(param = param, gradient = g)
      print(round(resdet,3))
      cat("--------------------------------------------\n")
    }
  }
  if (i >= iterlim) message = "maximum number of iterations reached"
  else message = "optimum reached"
  if (method == 'bfgs') H <- solve(Hm1)
  names(attr(x, 'gradient')) <- colnames(attr(x, 'gradi')) <- names(param)
  attr(x, "fixed") <- fixed
  est.stat = structure(list(elaps.time = NULL, nb.iter = i, eps = chi2,
    method = method, message = message), class = 'est.stat')
  result <- list(optimum = x,
                 coefficients = param,
                 est.stat = est.stat
                 )
  result
}


maxBFGSYC <- function(fn, grad=NULL, hess=NULL, start, print.level=0,
                  tol=1e-8, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  lambdatol=1e-6,
                  qrtol=1e-10,
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
   ## lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ## qrtol       - tolerance for qr decomposition
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
      "tol", "reltol", "gradtol", "steptol", "lambdatol", "qrtol",
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
       result <- maxNRCompute(fn=fn, grad=grad, hess=hess,
                              start=start,
                              print.level=print.level,
                              tol=tol, reltol=reltol,
                              gradtol=gradtol, steptol=steptol,
                              lambdatol=lambdatol,
                              qrtol=qrtol,
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
                        maxRoutine=maxNR,
                        constraints=constraints,
                        print.level=print.level,
                        tol=tol, reltol=reltol,
                        gradtol=gradtol, steptol=steptol,
                        lambdatol=lambdatol,
                        qrtol=qrtol,
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
