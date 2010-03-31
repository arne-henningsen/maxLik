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

   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim",
      "tol", "reltol", "temp", "tmax", "parscale" )
   checkFuncArgs( fn, argNames, "fn", "maxSANN" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxSANN" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxSANN" )
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
   hessian <- function(theta, ...) {
      ## just used for computing the final hessian, eventually using the supplied analytic information
      if(!is.null(hess)) {
         h <- as.matrix(hess(theta, ...))
      } else {
         if( is.null( grad ) ) {
            grad2 <- NULL
         } else {
            grad2 <- logLikGrad
         }
         h <- numericHessian( f = logLikFunc, grad = grad2, t0 = theta, ... )
      }
      rownames( h ) <- colnames( h ) <- names( start )
      return( h )
   }

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

   maximType <- "SANN maximisation"
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
      result <- optim(start, logLikFunc, control=control, method="SANN",
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
                          method="SANN", ...)
         resultConstraints <- list(type="constrOptim",
                                   barrier.value=result$barrier.value,
                                   outer.iterations=result$outer.iterations
                                   )
      }
      else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=fn, grad=grad, hess=hess,
                        start=start,
                        maxRoutine=maxSANN,
                        constraints=constraints,
                        print.level=print.level,
                        ...)
         return(result)
      }
      else {
         stop("maxSANN only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }
   result <- list(
                   maximum=result$value,
                   estimate=result$par,
                   gradient=logLikGrad(result$par, ...),
                   hessian=hessian(result$par, ...),
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

