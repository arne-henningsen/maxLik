## Return the #of parameters of model
nParam.maxim <- function(x, free=FALSE, ...) {
   if(!inherits(x, "maxim")) {
      stop("'nParam.maxim' called on non-'maxim' object")
   }
   return( .nParam( x, free=free ) )
}
