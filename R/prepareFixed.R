prepareFixed <- function( start, activePar, fixed ) {
   nParam <- length( start )
   ## establish the active parameters.
   if(!is.null(fixed)) {
      if(!is.null(activePar)) {
         if(!all(activePar)) {
            warning("Both 'activePar' and 'fixed' specified.  'activePar' ignored")
         }
      }
      activePar <- rep(TRUE, length(start))
      names(activePar) <- names(start)
      activePar[fixed] <- FALSE
   }
   else if(is.numeric(activePar)) {
      a <- rep(FALSE, nParam)
      a[activePar] <- TRUE
      activePar <- a
   }
   return( !activePar )
}