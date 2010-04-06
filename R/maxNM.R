maxNM <- function(fn, grad=NULL, hess=NULL,
                  start,
                  print.level=0,
                  iterlim=500,
                  constraints=NULL,
                  tol=1e-8, reltol=tol,
                  parscale=rep(1, length=length(start)),
                  alpha=1, beta=0.5, gamma=2,
                  ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   temp <- NULL
   tmax <- NULL

   method <- "Nelder-Mead"
   maxMethod <- "maxNM"

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
   environment( logLikHess ) <- environment()
   ## strip possible SUMT parameters and call the function thereafter
   environment( callWithoutSumt ) <- environment()
   maximType <- paste( method, "maximisation" )
   parscale <- rep(parscale, length.out=length(start))
   control <- list(trace=max(print.level, 0),
                   REPORT=1,
                   fnscale=-1,
                   reltol=reltol,
                   maxit=iterlim,
                   parscale=parscale,
                   alpha=alpha, beta=beta, gamma=gamma,
                   temp=temp, tmax=tmax )
   f1 <- callWithoutSumt( start, "logLikFunc", fnOrig = fn, gradOrig = grad, ... )
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
      result <- optim( par = start, fn = logLikFunc, control = control,
                      method = method, fnOrig = fn, gradOrig = grad, ... )
      resultConstraints <- NULL
   }
   else {
      ## linear equality and inequality constraints
                           # equality constraints: A %*% beta + B >= 0
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         ui <- constraints$ineqA
         ci <- -constraints$ineqB
         result <- constrOptim(theta=start, f=logLikFunc, grad=logLikGrad,
                          ui=ui, ci=ci, control=control,
                          method = method, fnOrig = fn, gradOrig = grad, ...)
         resultConstraints <- list(type="constrOptim",
                                   barrier.value=result$barrier.value,
                                   outer.iterations=result$outer.iterations
                                   )
      }
      else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=fn, grad=grad, hess=hess,
                        start=start,
                        maxRoutine = get( maxMethod ),
                        constraints=constraints,
                        print.level=print.level,
                        iterlim = iterlim,
                        tol = tol, reltol = reltol,
                        parscale = parscale,
                        alpha = alpha, beta= beta, gamma = gamma,
                        ...)
         return(result)
                           # this is already maxim object
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
   hessian <- logLikHess( result$par, fnOrig = fn, gradOrig = grad, ... )

   result <- list(
                  maximum=result$value,
                  estimate=result$par,
                  gradient=callWithoutSumt( result$par, "logLikGrad",
                     fnOrig = fn, gradOrig = grad, ... ),
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
