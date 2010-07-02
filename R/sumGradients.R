sumGradients <- function( gr, nParam ) {

   if( !is.null(dim(gr))) {
      gr <- colSums(gr)
   } else {
      ## ... or vector if only one parameter
      if(length(gr) > nParam ) {
         gr <- sum(gr)
      }
   }
   return( gr )
}