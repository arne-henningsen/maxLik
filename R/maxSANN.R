maxSANN <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=10000,
                    constraints = NULL,
                    tol=1e-8, reltol=tol,
                    temp=10, tmax=10, parscale=rep(1, length=length(start)),
                    ...) {
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
   ## We wrap the function in two versions:
   ## 1) sum over possible individual likelihoods
   ## 2) strip possible SUMT parameters and sum thereafter
   ## The former is for passing to the optimizer from withing sumt
   ## The latter is necessary for passing '...' to the function
   func <- function(theta, ...) {
      sum(fn(theta, ...))
   }
   funcS <- function(theta, ...) {
      ## this wrapper makes a) single-valued function (in case of BHHH
      ## vector-valued); and b) strips the SUMT extra arguments
      f <- match.call()
      f[names(formals(sumt))] <- NULL
      f[[1]] <- as.name("fn")
      names(f)[2] <- ""
      f1 <- eval(f, sys.frame(sys.parent()))
      sum(f1)
   }
   ## Needed for constrained optimization
   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         g <- grad(theta, ...)
         if(!is.null(dim(g))) {
            if(nrow(g) > 1) {
               g <- colSums( g )
            }
         }
         names( g ) <- names( start )
         return( g )
      }
      g <- numericGradient(fn, theta, ...)
      if(!is.null(dim(g))) {
         return(colSums(g))
      } else {
         return(g)
      }
   }
   hessian <- function(theta, ...) {
      ## just used for computing the final hessian, eventually using the supplied analytic information
      if(!is.null(hess)) {
         h <- as.matrix(hess(theta, ...))
      } else {
         h <- numericHessian(func, gradient, theta, ...)
      }
      rownames( h ) <- colnames( h ) <- names( start )
      return( h )
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
                   tmax=tmax)
   f1 <- funcS(start, ...)
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
      result <- optim(start, func, control=control, method="SANN",
                      hessian=FALSE, ...)
      resultConstraints <- NULL
   }
   else {
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         ui <- constraints$ineqA
         ci <- -constraints$ineqB
         result <- constrOptim(theta=start, f=func, grad=gradient,
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
                   gradient=gradient(result$par, ...),
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

