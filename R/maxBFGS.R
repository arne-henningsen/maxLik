maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=200,
                    constraints=NULL,
                    tol=1e-8, reltol=tol,
                    parscale=rep(1, length=length(start)),
                    ## sumt parameters
                    ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ...           further arguments to fn() and grad()

   alpha <- NULL
   beta <- NULL
   gamma <- NULL
   temp <- NULL
   tmax <- NULL

   method <- "BFGS"
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
   environment( logLikGrad ) <- environment()
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
   f1 <- callWithoutSumt( start, "logLikFunc", func = fn, ...)
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
   G1 <- callWithoutSumt( start, "logLikGrad", func = fn, ...)
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != length(start)) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", length(start), ")" )
   }
   ## A note about return value:
   ## We can the return from 'optim' in a object of class 'maxim'.
   ## However, as 'sumt' already returns such an object, we return the
   ## result of 'sumt' directly, without the canning
   if(is.null(constraints)) {
       result <- optim( par = start, fn = logLikFunc, control = control,
                      method = method, gr = logLikGrad, func = fn, ... )
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
                          method = method, func = fn, ...)
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
   hessian <- logLikHess( result$par, func = fn, ... )

   result <- list(
                   maximum=result$value,
                   estimate=result$par,
                   gradient=logLikGrad( theta = result$par, func = fn, ... ),
                   hessian=hessian,
                   code=result$convergence,
                   message=paste(message(result$convergence), result$message),
                   last.step=NULL,
                   activePar = rep( TRUE, length ( result$par ) ),
                   iterations=result$counts,
                   type=maximType,
                  constraints=resultConstraints
                  )
   class(result) <- "maxim"
   invisible(result)
}
