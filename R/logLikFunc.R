## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, sumObs = TRUE, ...) {

   # Arguments "gradOrig" and "hessOrig" are just for compatibility with
   #    logLikGrad() and logLikHess()

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   result <- fnOrig( theta, ... )

   if( sumObs ) {
      result <- sumKeepAttr( result )
      g <- attributes( result )$gradient
      if( !is.null( g ) ) {
         g <- sumGradients( g, length( theta ) )
         names( g ) <- names( theta )
         if( !is.null( fixed ) ) {
            g <- g[ !fixed ]
         }
         attributes( result )$gradient <- g
      }
   }

   return( result )
}
