## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, gradOrig, hessOrig, ...) {
   if(is.null(gradOrig)) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig, ...)
   } else {
      g <- gradOrig(theta, ...)
   }
   if(!is.null(dim(g))) {
      g <- colSums(g)
   }
   names( g ) <- names( theta )
   return( g )
}
