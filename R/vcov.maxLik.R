## maxLik
vcov.maxLik <- function(object, eigentol=1e-12, ...) {
   ## if exists $varcovar, take it
   if(!is.null(object$varcovar))
       return(object$varcovar)
   ## otherwise invert hessian
   activePar <- activePar(object)
   hess <- hessian(object)[activePar, activePar] 
   hessev <- abs(eigen(hess, symmetric=TRUE, only.values=TRUE)$values)
   if(min(hessev) > (eigentol*max(hessev))) {
      varcovar <- matrix(0, nParam(object), nParam(object))
      varcovar[activePar,activePar] <- solve(-hessian(object)[activePar,activePar])
   }
   else
       varcovar <- NULL
   varcovar
}
