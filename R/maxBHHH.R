maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods
   func <- function(theta, ...) {
      ## we wrap the likelihood function here, in order to save gradient
      ## value when it is supplied as attribute
      a <- fn(theta, ...)
      if(!is.null(grad <- attr(a, "gradient")))
          assign("gradVal", grad, inherits=TRUE)
      attr(a, "hessian") <- NULL
                           # this is to ensure we do BHHH even if the
                           # likelihood function supplies analytic Hessian.
                           # Use maxNR if you have Hessian and don't
                           # want to do BHHH!
      a
   }
   ##
   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim" )
   checkFuncArgs( fn, argNames, "fn", "maxBHHH" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBHHH" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBHHH" )
   }
   ## Save the value of gradient and use it later for hessian
   gradVal <- NULL


   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         ## Gradient supplied by the user.
         g <- grad(theta, ...)
         if( length( theta ) == 1 ) {
            g <- matrix( g, ncol = 1 )
         }
      } else {
         ## fall back to the numeric gradient
         g <- numericGradient(fn, theta, ...)
      }
      checkBhhhGrad( g = g, theta = theta, analytic = !is.null( grad ) )
      assign("gradVal", g, inherits=TRUE)
      return( g )
   }
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }
   a <- maxNR(func, grad=gradient, hess=hess, start=start,
              iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
