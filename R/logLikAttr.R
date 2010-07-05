logLikAttr <- function(theta, fnOrig, gradOrig, hessOrig, fixed,
         sumObs = FALSE, ...) {

# this function returns the log-likelihood value with gradient and Hessian as
# attributes. If the log-likelihood function provided by the user does not add
# these attributes, this functions uses the functions provided by the user
# as arguments "grad" and "hess" or (if they are not provided) uses the
# finite-difference method to obtain the gradient and Hessian

         # large initial indentation to be able to diff to previous version
         # that was defined in maxNR() / maxNR.R.

         ## number of parameters
         nParam <- length( theta )

         ## value of log-likelihood function
         f <- fnOrig(theta, ...)
         if( sumObs ) {
            fAttr <- attributes(f)
            f <- sum(f)
            mostattributes(f) <- fAttr
         }

         ## gradient of log-likelihood function
         gr <- attr( f, "gradient" )
         if( is.null( gr ) ) {
            if( !is.null( gradOrig ) ) {
               gr <- gradOrig(theta, ...)
            } else {
               gr <- numericGradient(f = fnOrig, t0 = theta,
                                    fixed=fixed, ...)
            }
         }
         ## Now check if the gradient is vector or matrix...
         if(!sumObs) {
            if(observationGradient(gr, length(theta))) {
               gr <- as.matrix(gr)
            }
         } else {
            ## We need just summed gradient
            gr <- sumGradients( gr, nParam )
         }
         ## Set gradients of fixed parameters to zero so that they are always zero
         ## (no matter if they are analytical or finite-difference gradients)
         if( is.null( dim( gr ) ) ) {
            gr[ fixed ] <- NA
         } else {
            gr[ , fixed ] <- NA
         }

         ## Hessian of log-likelihood function
         h <- attr( f, "hessian" )
         if( is.null( h ) ) {
            if(!is.null(hessOrig)) {
               h <- as.matrix(hessOrig(theta, ...))
            } else {
               llFunc <- function( theta, ... ) {
                  return( sum( fnOrig( theta, ... ) ) )
               }
               if( !is.null( attr( f, "gradient" ) ) ) {
                  gradFunc <- function( theta, ... ) {
                     return( sumGradients( attr( fnOrig( theta, ... ), "gradient" ),
                        nParam ) )
                  }
               } else if( !is.null( gradOrig ) ) {
                  gradFunc <- function( theta, ... ) {
                     return( sumGradients( gradOrig( theta, ... ), nParam ) )
                  }
               } else {
                  gradFunc <- NULL
               }
               h <- numericHessian( f = llFunc, grad = gradFunc, t0 = theta,
                                 fixed=fixed, ...)
            }
         }
         if((dim(h)[1] != nParam) | (dim(h)[2] != nParam)) {
            stop("Wrong hessian dimension.  Needed ", nParam, "x", nParam,
               " but supplied ", dim(h)[1], "x", dim(h)[2])
         }
         ## Set elements of the Hessian corresponding to the fixed parameters
         ## to zero so that they are always zero (no matter if they are
         ## calculated analytical or by the finite-difference method)
         h[ fixed, ] <- NA
         h[ , fixed ] <- NA

         attr( f, "gradient" ) <- gr
         attr( f, "hessian" ) <- h
         return( f )
}
