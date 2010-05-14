## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, sumObs = TRUE, ...) {

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   if(is.null(gradOrig)) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig,
         sumObs = sumObs, ...)
   } else {
      g <- gradOrig(theta, ...)
   }
   if( sumObs ) {
      if(!is.null(dim(g)) ) {
         g <- colSums(g)
      } else if( length( theta ) == 1 && length( g ) > 1 ) {
         g <- sum( g )
      }
   }
   if(is.null(dim(g))) {
      names( g ) <- names( theta )
      if( !is.null( fixed ) ) {
         g <- g[ !fixed ]
      }
   } else {
      colnames( g ) <- names( theta )
      if( !is.null( fixed ) ) {
         g <- g[ , !fixed ]
      }
   }

   return( g )
}
