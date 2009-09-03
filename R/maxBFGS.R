maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=200,
                    constraints,
                    tol=1e-8, reltol=tol,
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
   func <- function(theta, ...) {
      sum(fn(theta, ...))
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
      g <- numericGradient(fn, theta, ...)
      if(!is.null(dim(g))) {
         return(colSums(g))
      } else {
         return(g)
      }
   }
   type <- "BFGS maximisation"
   control <- list(trace=print.level,
                    REPORT=1,
                    fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim)
   if(missing(constraints))
       a <- optim(start, func, gr=gradient, control=control, method="BFGS",
                  ...)
   else {
      ui <- constraints$ui
      ci <- constraints$ci
      a <- constrOptim(theta=start, f=func, grad=gradient, ui=ui, ci=ci, control=control,
                       method="BFGS", ...)
   }
   if(!is.null(hess)) {
       hessian <- hess(a$par)
   } else {
      if( is.null( grad ) ) {
         grad2 <- NULL
      } else {
         grad2 <- gradient
      }
      hessian <- numericHessian( f = func, grad = grad2, t0=a$par, ... )
   }
   rownames( hessian ) <- colnames( hessian ) <- names( a$par )
   result <- list(
                   maximum=a$value,
                   estimate=a$par,
                   gradient=gradient( theta = a$par, ... ),
                   hessian=hessian,
                   code=a$convergence,
                   message=paste(message(a$convergence), a$message),
                   last.step=NULL,
                   activePar = rep( TRUE, length ( a$par ) ),
                   iterations=a$counts,
                   type=type)
   class(result) <- "maxim"
   invisible(result)
}

