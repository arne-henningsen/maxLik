maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=200,
                    constraints=NULL,
                    tol=1e-8, reltol=tol,
                    ## sumt parameters
                    ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ...           further arguments to fn() and grad()
   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim",
      "constraints", "tol", "reltol" )
   checkFuncArgs( fn, argNames, "fn", "maxBFGS" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBFGS" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBFGS" )
   }

   message <- function(c) {
      switch(as.character(c),
               "0" = "successful convergence",
               "10" = "degeneracy in Nelder-Mead simplex"
               )
   }
   ##
   ## sum over possible individual likelihoods or gradients
   environment( logLikFunc ) <- environment()
   environment( logLikGrad ) <- environment()
   ## strip possible SUMT parameters and call the function thereafter
   environment( callWithoutSumt ) <- environment()
   maximType <- "BFGS maximisation"
   control <- list(trace=print.level,
                    REPORT=1,
                    fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim)
   f1 <- callWithoutSumt( start, "logLikFunc", ...)
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
   G1 <- callWithoutSumt( start, "logLikGrad", ...)
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
       result <- optim(start, logLikFunc, gr=logLikGrad, control=control,
                       method="BFGS",
                       ...)
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
                          method="BFGS", ...)
         resultConstraints <- list(type="constrOptim",
                                   barrier.value=result$barrier.value,
                                   outer.iterations=result$outer.iterations
                                   )
      }
      else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         result <- sumt(fn=fn, grad=grad, hess=hess,
                        start=start,
                        maxRoutine=maxBFGS,
                        constraints=constraints,
                        print.level=print.level,
                        ...)
         return(result)
                           # this is already maxim object
      }
      else {
         stop("maxBFGS only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }
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
                   gradient=logLikGrad( theta = result$par, ... ),
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
