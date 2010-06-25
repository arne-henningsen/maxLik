## maxLik
vcov.maxLik <- function(object, eigentol=1e-12, ...) {
   ## if exists $varcovar, take it
   if(!is.null(object$varcovar))
       return(object$varcovar)
   ## otherwise invert hessian
   activePar <- activePar(object)
   if(!is.null(hess <- hessian(object))) {
      hess <- hessian(object)[activePar, activePar,drop=FALSE] 
      hessev <- abs(eigen(hess, symmetric=TRUE, only.values=TRUE)$values)
      varcovar <- matrix(0, nParam(object), nParam(object))
                           # this makes the fixed parameters to 0
      rownames( varcovar ) <- colnames(varcovar ) <- names(coef(object))
      if(min(hessev) > (eigentol*max(hessev))) {
      ## If hessian is not singular, fill in the free parameter values
         varcovar[activePar,activePar] <- solve(-hessian(object)[activePar,activePar])
      }
      else {
      ## If singular, the free parameter values will be Inf
         varcovar[activePar,activePar] <- Inf
      }
      return(varcovar)
   }
   else
       return(NULL)
}
