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
      g <- numericGradient(func, theta, ...)
      if(!is.null(dim(g))) {
         return(colSums(g))
      } else {
         return(g)
      }
   }
   gradientS <- function(theta, ...) {
      if(!is.null(grad)) {
         g <- match.call()
         g[names(formals(sumt))] <- NULL
         g[[1]] <- as.name("grad")
         names(g)[2] <- ""
         g <- eval(g, sys.frame(sys.parent()))
         if(!is.null(dim(g))) {
            if(nrow(g) > 1) {
               g <- colSums( g )
            }
         }
         names( g ) <- names( start )
         return( g )
      }
      g <- match.call()
      g[names(formals(sumt))] <- NULL
      g[[1]] <- as.name("numericGradient")
      names(g)[2] <- "t0"
      g$f <- func
      g <- eval(g, sys.frame(sys.parent()))
      if(!is.null(dim(g))) {
         return(colSums(g))
      } else {
         return(g)
      }
   }
   maximType <- "BFGS maximisation"
   control <- list(trace=print.level,
                    REPORT=1,
                    fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim)
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
   G1 <- gradientS(start, ...)
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
       result <- optim(start, func, gr=gradient, control=control,
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
         result <- constrOptim(theta=start, f=func, grad=gradient,
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
         grad2 <- gradient
      }
      hessian <- numericHessian( f = func, grad = grad2, t0=result$par, ... )
   }
   rownames( hessian ) <- colnames( hessian ) <- names( result$par )
   result <- list(
                   maximum=result$value,
                   estimate=result$par,
                   gradient=gradient( theta = result$par, ... ),
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
