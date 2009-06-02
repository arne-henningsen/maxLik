maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods

   ## checking the log likelihodd function
   testVal <- fn( start, ... )
   if( length( testVal ) < length( start ) ) {
      stop( "if the gradients (argument 'grad') are not provided by the user,",
         " the BHHH method requires that the log-likelihood function",
         " (argument 'fn') returns a numeric vector with at least as many",
         " elements as the number of parameters (", length( start ), "),",
         " where each element must be the log-likelihood value corresponding",
         " to an individual observation" )
   }
   rm( testVal )

   ## checking the gradient function
   if( !is.null( grad ) ) {
      testGrad <- grad( start, ... )
      if( nrow( matrix( testGrad ) ) < length( start ) ) {
         stop( "the BHHH method requires that the gradient function",
            " (argument 'grad') returns a numeric matrix with at least as many",
            " rows as the number of parameters (", length( start ), "),",
            " where each row must correspond to the gradients",
            " of the log-likelihood function of an individual observation" )
      }
      rm( testGrad )
   }

   ## Save the value of gradient and use it later for hessian
   gradVal <- NULL

   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         g <- grad(theta, ...)
      } else {
         g <- numericGradient(fn, theta, ...)
      }
      ## Ensure g is suitable for information equality approximation
      if(is.null(dim(g)))
          g <- matrix(g)
      if(!(dim(g)[1] > length(theta) & dim(g)[2] == length(theta))) {
         stop(paste("Gradient matrix must have at least as many rows and exactly as many columns as the number of parameters.\n",
                    "Currently", length(theta), "parameters but the gradient is", dim(g)[1], "x", dim(g)[2]))
      }
      ##
      assign("gradVal", g, inherits=TRUE)
      return( g )
   }
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }
   a <- maxNR(fn, grad=gradient, hess=hess, start=start, iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
