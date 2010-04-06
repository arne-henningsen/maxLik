## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, ...) {
   if(is.null(grad)) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig, ...)
   } else {
      g <- grad(theta, ...)
   }
   if(!is.null(dim(g))) {
      g <- colSums(g)
   }
   names( g ) <- names( theta )
   return( g )
}
