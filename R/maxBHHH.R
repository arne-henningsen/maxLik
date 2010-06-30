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

   ## extract attributes "gradient" and "hessian"
   f1 <- callWithoutArgs( theta = start, fName = "fn", 
      args = c( names( formals( maxNR ) ), names( formals(sumt) ) ), ... )
   attGradient <- attributes( f1 )$gradient
   attHessian  <- attributes( f1 )$hessian

   ## initialize variable for saving the value of gradient 
   ## so that it can be used later for calculating the Hessian
   gradVal <- NULL

   ## wrapper function for obtaining the log-likelihood value
   func <- function(theta, ...) {
      ## we wrap the likelihood function here, in order to save gradient
      ## value when it is supplied as attribute
      a <- fn(theta, ...)
      if(!is.null(grad <- attr(a, "gradient"))) {
         checkBhhhGrad( g = grad, theta = theta, analytic = TRUE )
          assign("gradVal", grad, inherits=TRUE)
      }
      attr(a, "hessian") <- NULL
                           # this is to ensure we do BHHH even if the
                           # likelihood function supplies analytic Hessian.
                           # Use maxNR if you have Hessian and don't
                           # want to do BHHH!
      a
   }

   ## wrapper function for obtaining the gradients
   if( is.null( attGradient ) ) {
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
   } else {
      gradient <- NULL
   }

   ## wrapper function for obtaining the Hessian
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }

   ## using the Newton-Raphson algorithm with BHHH method for Hessian
   a <- maxNR(func, grad=gradient, hess=hess, start=start,
              iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
