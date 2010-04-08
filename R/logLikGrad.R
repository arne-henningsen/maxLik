## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, ...) {

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   if(is.null(gradOrig)) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig, ...)
   } else {
      g <- gradOrig(theta, ...)
   }
   if(!is.null(dim(g))) {
      g <- colSums(g)
   }
   names( g ) <- names( theta )

   if( !is.null( fixed ) ) {
      g <- g[ !fixed ]
   }

   return( g )
}
