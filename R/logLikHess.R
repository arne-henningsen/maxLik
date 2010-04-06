logLikHess <- function( theta, fnOrig, gradOrig, ... ) {
   if(!is.null(hess)) {
       hessian <- hess( theta, ... )
   } else {
      if( is.null( gradOrig ) ) {
         grad2 <- NULL
      } else {
         grad2 <- logLikGrad
      }
      hessian <- numericHessian( f = logLikFunc, grad = grad2, t0 = theta, 
         fnOrig = fnOrig, gradOrig = gradOrig, ... )
   }
   rownames( hessian ) <- colnames( hessian ) <- names( theta )
   return( hessian )
}
