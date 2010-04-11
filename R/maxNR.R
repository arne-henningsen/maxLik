maxNR <- function(fn, grad=NULL, hess=NULL, start, print.level=0,
                  tol=1e-8, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  lambdatol=1e-6,
                  qrtol=1e-10,
                  iterlim=150,
                  constraints=NULL,
                  activePar=NULL,
                  fixed=NULL,
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
