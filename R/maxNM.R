maxNM <- function(fn, grad=NULL, hess=NULL,
                  start,
                  print.level=0,
                  iterlim=500,
                  constraints,
                  tol=1e-8, reltol=tol,
                  parscale=rep(1, length=length(start)),
                  alpha=1, beta=0.5, gamma=2,
                  ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim",
      "constraints", "tol", "reltol", "parscale", "alpha", "beta", "gamma" )
   checkFuncArgs( fn, argNames, "fn", "maxNM" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxNM" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxNM" )
   }

   message <- function(c) {
      switch(as.character(c),
             "0" = "successful convergence",
             "10" = "degeneracy in Nelder-Mead simplex",
             "51" = "warning from the 'L-BFGS-B' method; see the corresponding component 'message' for details",
             "52" = "error from the 'L-BFGS-B' method; see the corresponding component 'message' for details"
             )
   }
   func <- function(theta, ...) {
      sum(fn(theta, ...))
   }
   gradient <- function(theta, ...) {
      ## this gradient function is used only for returning the gradient at the optimum
      if(!is.null(grad)) {
         g <- grad(theta, ...)
         if(!is.null(dim(g))) {
            if(ncol(g) > 1) {
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
         h <- numericHessian(fn, gradient, theta, ...)
      }
      rownames( h ) <- colnames( h ) <- names( start )
      return( h )
   }
   type <- "Nelder-Mead maximisation"
   parscale <- rep(parscale, length.out=length(start))
   control <- list(trace=print.level,
                   REPORT=1,
                   fnscale=-1,
                   reltol=reltol,
                   maxit=iterlim,
                   parscale=parscale,
                   alpha=alpha, beta=beta, gamma=gamma
                   )
   if(missing(constraints))
       a <- optim(start, func, control=control, method="Nelder-Mead", hessian=FALSE, ...)
   else {
      ui <- constraints$ui
      ci <- constraints$ci
      a <- constrOptim(theta=start, f=func, ui=ui, ci=ci, control=control,
                       method="Nelder-Mead", ...)
   }
   result <- list(
                  maximum=a$value,
                  estimate=a$par,
                  gradient=gradient(a$par),
                  hessian=hessian(a$par),
                  code=a$convergence,
                  message=paste(message(a$convergence), a$message),
                  last.step=NULL,
                  activePar = rep( TRUE, length ( a$par ) ),
                  iterations=a$counts[1],
                  type=type)
   class(result) <- "maxim"
   invisible(result)
}
