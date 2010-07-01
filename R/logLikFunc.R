## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, sumObs = TRUE, ...) {

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   result <- fnOrig( theta, ... )

   if( sumObs ) {
      resultAttr <- attributes( result )
      result <- sum( result )
      mostattributes( result ) <- attributes( result )
      g <- attributes( result )$gradient
      if( !is.null( g ) ) {
         if( !is.null( dim( g ) ) ) {
            g <- colSums(g)
         } else if( length( theta ) == 1 && length( g ) > 1 ) {
            g <- sum( g )
         }
         names( g ) <- names( theta )
         if( !is.null( fixed ) ) {
            g <- g[ !fixed ]
         }
         attributes( result )$gradient <- g
      }
   }

   return( result )
}
