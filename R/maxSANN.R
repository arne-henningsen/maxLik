maxSANN <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=10000,
                    constraints = NULL,
                    tol=1e-8, reltol=tol,
                    temp=10, tmax=10, parscale=rep(1, length=length(start)),
                    random.seed = 123, ... ) {
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   method <- "SANN"
   maxMethod <- paste( "max", method, sep = "" )

   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim",
      "constraints", "tol", "reltol", "parscale", "alpha", "beta", "gamma",
      "temp", "tmax" )
   checkFuncArgs( fn, argNames, "fn", maxMethod )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", maxMethod )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", maxMethod )
   }

   message <- function(c) {
      switch(as.character(c),
             "0" = "successful convergence",
             "10" = "degeneracy in Nelder-Mead simplex",
             "51" = "warning from the 'L-BFGS-B' method; see the corresponding component 'message' for details",
             "52" = "error from the 'L-BFGS-B' method; see the corresponding component 'message' for details"
             )
   }
   ##
   ## sum over possible individual likelihoods or gradients
   environment( logLikFunc ) <- environment()
   environment( logLikGrad ) <- environment()
   ## strip possible SUMT parameters and call the function thereafter
   environment( callWithoutSumt ) <- environment()

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by 'optim( method="SANN" )')
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   maximType <- paste( method, "maximisation" )
   parscale <- rep(parscale, length.out=length(start))
   control <- list(trace=print.level,
                    REPORT=1,
                   fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim,
                   parscale=parscale,
                   temp=temp,
                   tmax=tmax )
   f1 <- callWithoutSumt( start, "logLikFunc", ... )
   if(is.na( f1)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maximType)
      class(result) <- "maxim"
      return(result)
   }
   if(print.level > 2) {
      cat("Initial function value:", f1, "\n")
   }
   ## A note about return value:
   ## We can the return from 'optim' in a object of class 'maxim'.
   ## However, as 'sumt' already returns such an object, we return the
   ## result of 'sumt' directly, without the canning
   if(is.null(constraints)) {
      result <- optim(start, logLikFunc, control=control, method = method,
                      hessian=FALSE, ...)
      resultConstraints <- NULL
   }
   else {
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         ui <- constraints$ineqA
         ci <- -constraints$ineqB
         result <- constrOptim(theta=start, f=logLikFunc,
                               # Note that gradient has different meaning for SANN!
                          ui=ui, ci=ci, control=control,
                          method = method, ...)
         resultConstraints <- list(type="constrOptim",
                                   barrier.value=result$barrier.value,
                                   outer.iterations=result$outer.iterations
                                   )
      }
      else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=fn, grad=grad, hess=hess,
                        start=start,
                        maxRoutine=get( maxMethod ),
                        constraints=constraints,
                        print.level=print.level,
                        ...)
         return(result)
      }
      else {
         stop( maxMethod, " only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }

   # calculate (final) Hessian
   if(!is.null(hess)) {
       hessian <- hess(result$par)
   } else {
      if( is.null( grad ) ) {
         grad2 <- NULL
      } else {
         grad2 <- logLikGrad
      }
      hessian <- numericHessian( f = logLikFunc, grad = grad2, t0=result$par, ... )
   }
   rownames( hessian ) <- colnames( hessian ) <- names( result$par )

   result <- list(
                   maximum=result$value,
                   estimate=result$par,
                   gradient=logLikGrad(result$par, ...),
                   hessian=hessian,
                   code=result$convergence,
                   message=paste(message(result$convergence), result$message),
                   last.step=NULL,
                   activePar = rep( TRUE, length ( result$par ) ),
                   iterations=result$counts[1],
                   type=maximType,
                  constraints=resultConstraints
                  )
   class(result) <- "maxim"
   invisible(result)
}

