logLikHess <- function( theta, fnOrig, gradOrig, hessOrig,
      start = NULL, fixed = NULL, gradAttr = NULL, hessAttr = NULL, ... ) {
   ## Calculate the Hessian of the function, either by analytic or numeric method

   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)
   
   if(!is.null(hessOrig)) {
       hessian <- as.matrix(hessOrig( theta, ... ))
   } else {
      if( is.null( hessAttr ) || hessAttr || is.null( gradAttr ) ) {
         llVal <- fnOrig( theta, ... )
         gradient <- attr( llVal, "gradient" )
         hessian <- attr( llVal, "hessian" )
         gradAttr <- !is.null( gradient )
         hessAttr <- !is.null( hessian )
      }
      if( !hessAttr ) {
         if( !is.null( gradOrig ) ) {
            grad2 <- logLikGrad
         } else if( gradAttr ) {
            grad2 <- function( theta, fnOrig = NULL, gradOrig = NULL, ... ) {
               gradient <- attr( fnOrig( theta, ... ), "gradient" )
               gradient <- sumGradients( gradient, length( theta ) )
               return( gradient )
            }
         } else {
            grad2 <- NULL
         }
         hessian <- numericHessian( f = logLikFunc, grad = grad2, t0 = theta,
            fnOrig = fnOrig, gradOrig = gradOrig, ... )
      }
   }
   rownames( hessian ) <- colnames( hessian ) <- names( theta )

   if( !is.null( fixed ) ) {
      hessian <- hessian[ !fixed, !fixed, drop = FALSE ]
   }

   return( hessian )
}
