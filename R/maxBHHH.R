maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods

   ## check if arguments of user-provided functions have reserved names
   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim" )
   checkFuncArgs( fn, argNames, "fn", "maxBHHH" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBHHH" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBHHH" )
   }

   ## wrapper function for obtaining the log-likelihood value, gradient, and Hessian
   func <- function(theta, ...) {
      a <- fn(theta, ...)
      if(!is.null(g <- attr(a, "gradient"))) {
         ## gradient when it is supplied as attribute
         checkBhhhGrad( g = g, theta = theta, analytic = TRUE )
      } else if(!is.null(grad)) {
         ## gradient supplied as function "grad"
         g <- grad(theta, ...)
         if( length( theta ) == 1 ) {
            g <- matrix( g, ncol = 1 )
         }
         checkBhhhGrad( g = g, theta = theta, analytic = TRUE )
      } else {
         ## fall back to the numeric gradient
         g <- numericGradient(fn, theta, ...)
         checkBhhhGrad( g = g, theta = theta, analytic = FALSE )
      }
      attr( a, "gradient" ) <- g
      attr( a, "hessian" ) <- -t( g ) %*% g

      a
   }

   ## using the Newton-Raphson algorithm with BHHH method for Hessian
   a <- maxNR(func, start=start,
              iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
