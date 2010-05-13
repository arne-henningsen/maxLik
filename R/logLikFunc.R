## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, sumObs = TRUE, ...) {

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   result <- fnOrig( theta, ... )

   if( sumObs ) {
      result <- sum( result )
   }

   return( result )
}
