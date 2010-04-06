## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, func, ...) {
   if(is.null(grad)) {
      g <- numericGradient(logLikFunc, theta, func = func, ...)
   } else {
      g <- grad(theta, ...)
   }
   if(!is.null(dim(g))) {
      g <- colSums(g)
   }
   names( g ) <- names( start )
   return( g )
}
