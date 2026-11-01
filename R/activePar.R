## activePar: returns parameters which are free under maximisation (not fixed as constants)

activePar <- function(x, ...)
    UseMethod("activePar")

activePar.default <- function(x, ...) {
   result <- .activePar( x )
   if( is.null( result ) ) {
      ## fallback for arbitrary objects
      result <- rep( TRUE, length( coef( x ) ) )
   }
   return( result )
}
