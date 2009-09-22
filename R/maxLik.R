maxLik <- function(logLik, grad=NULL, hess=NULL, start,
                   method="Newton-Raphson",
                   constraints,
                   ...) {
   ## Maximum Likelihood estimation.
   ##
   ## Newton-Raphson maximisation
   ## Parameters:
   ## logLik     log-likelihood function.  First argument must be the vector of parameters.
   ## grad       gradient of log-likelihood.  If NULL, numeric gradient is used.  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(nObs, nParam).  In this case the rows are simply
   ##                 summed (useful for maxBHHH).
   ## hess       Hessian function (numeric used if NULL)
   ## start      initial vector of parameters (eventually w/names)
   ## method     maximisation method (Newton-Raphson)
   ## constraints  constrained optimization: a list (see below)
   ## ...        additional arguments for the maximisation routine
   ##
   ## RESULTS:
   ## list of class c("maxLik", "maxim").  This is in fact equal to class "maxim", just the
   ## methods are different.
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success, depends on the optimization
   ##             method 
   ## message     character message describing the code
   ## type        character, type of optimization
   ##
   ##             there may be more components, depending on the choice of
   ##             the algorith.
   ##             
   argNames <-  c( "logLik", "grad", "hess", "start", "method",
                  "constraints" )
   checkFuncArgs( logLik, argNames, "logLik", "maxLik" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxLik" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxLik" )
   }
   maxRoutine <- switch(method,
                        "Newton-Raphson" =,
                        "newton-raphson" =,
                        "NR" =,
                        "nr" = maxNR,
                        "BFGS" =,
                        "bfgs" = maxBFGS,
                        "BHHH" =,
                        "bhhh" = maxBHHH,
                        "Nelder-Mead" =,
                        "NM" =,
                        "nm" = maxNM,
                        "SANN" =,
                        "sann" = maxSANN,
                        stop( "Maxlik: unknown maximisation method ", method )
                        )
   ## Constrained optimization.  We can two possibilities:
   ## * linear equality constraints
   ## * linear inequality constraints
   ##
   if(!missing(constraints)) {
      if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=logLik, grad=grad, hess=hess,
                        start=start,
                        maxRoutine=maxRoutine,
                        constraints=constraints, ...)
      }
      else if(identical(names(constraints), c("ineqA", "ineqB"))) {
                           # inequality constraints A %*% beta + B > 0
      }
      else stop("'constraints' must be a list with two components:",
                "either 'eqA' and 'eqB' or 'ineqA' and 'ineqB' for",
                "equality and inequality constraints respectively")
   }
   else
                           # unconstrained optimization
       result <- maxRoutine(fn=logLik, grad=grad, hess=hess, start=start,
                            ...)
   class(result) <- c("maxLik", class(result))
   result
}
