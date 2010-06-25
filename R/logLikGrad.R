## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, gradOrig, hessOrig,
                       start = NULL, fixed = NULL, sumObs = TRUE,
                       suppliedValue=NULL,
                       ...) {
   ## suppliedValue   a call to the function may provide the (pre-calculated) value of gradient.  If non-NULL, this is used.
   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)
   if(!is.null(suppliedValue)) {
      g <- suppliedValue
   }
   else if(is.null(gradOrig)) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig,
         sumObs = sumObs, ...)
   } else {
      g <- gradOrig(theta, ...)
   }
   if( sumObs ) {
      ## We were requested a single (summed) gradient.  Return a vector
      if(!is.null(dim(g)) ) {
         g <- colSums(g)
      } else if( length( theta ) == 1 && length( g ) > 1 ) {
         g <- sum( g )
      }
      names( g ) <- names( theta )
      if( !is.null( fixed ) ) {
         g <- g[ !fixed ]
      }
   }
   else {
      ## we were requested individual gradients (if possible).  Ensure g is a matrix
      if(observationGradient(g, length(theta))) {
         ## it was indeed by observations
         g <- as.matrix(g)
         colnames( g ) <- names( theta )
         if( !is.null( fixed ) ) {
            g <- g[ , !fixed ]
         }
      }
      else {
         ## it wasn't
         g <- drop(g)
         names(g) <- names(theta)
         if( !is.null( fixed ) ) {
            g <- g[ !fixed ]
         }
      }         
   }
   return( g )
}
