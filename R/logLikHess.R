logLikHess <- function( theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, ... ) {
   ## Calculate the Hessian of the function, either by analytic or numeric method

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)
   
   if(!is.null(hessOrig)) {
       hessian <- as.matrix(hessOrig( theta, ... ))
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

   if( !is.null( fixed ) ) {
      hessian <- hessian[ !fixed, !fixed, drop = FALSE ]
   }

   return( hessian )
}
