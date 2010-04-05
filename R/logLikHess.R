logLikHess <- function( theta, ... ) {
   if(!is.null(hess)) {
       hessian <- hess( theta )
   } else {
      if( is.null( grad ) ) {
         grad2 <- NULL
      } else {
         grad2 <- logLikGrad
      }
      hessian <- numericHessian( f = logLikFunc, grad = grad2, t0 = theta, ... )
   }
   rownames( hessian ) <- colnames( hessian ) <- names( theta )
   return( hessian )
}
