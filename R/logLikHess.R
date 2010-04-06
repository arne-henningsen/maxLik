logLikHess <- function( theta, fnOrig, gradOrig, hessOrig, ... ) {
   if(!is.null(hessOrig)) {
       hessian <- hessOrig( theta, ... )
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
