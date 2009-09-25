maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods

   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim" )
   checkFuncArgs( fn, argNames, "fn", "maxBHHH" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBHHH" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBHHH" )
   }

   ## Save the value of gradient and use it later for hessian
   gradVal <- NULL

   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         ## Gradient supplied by the user.
         g <- grad(theta, ...)
         if( length( theta ) == 1 ) {
            g <- matrix( g, ncol = 1 )
         }
         ## Check whether the gradient has enough rows (about enough
         ## observations in data)
         if( !is.matrix( g ) ) {
            stop("Gradient is not a matrix.\n",
                 "the BHHH method requires that the gradient function\n",
                 "(argument 'grad') returns a numeric matrix,\n",
                 "where each row must correspond to the gradient(s)\n",
                 "of the log-likelihood function at an individual\n",
                 "(independent) observation and each column must\n",
                 "correspond to a parameter" )
         } else if( nrow( g ) < length( theta ) ) {
            stop( "the matrix returned by the gradient function",
               " (argument 'grad') must have at least as many",
               " rows as the number of parameters (", length( theta ), "),",
               " where each row must correspond to the gradients",
               " of the log-likelihood function of an individual",
               " (independent) observation:\n",
               " currently, there are (is) ", length( theta ), " parameter(s)",
               " but the gradient matrix has only ", nrow( g ), " row(s)" )
         } else if( ncol( g ) != length( theta ) ) {
            stop( "the matrix returned by the gradient function",
               " (argument 'grad') must have exactly as many columns",
               " as the number of parameters:\n",
               " currently, there are (is) ", length( theta ), " parameter(s)",
               " but the gradient matrix has ", ncol( g ), " columns" )
         }
      } else {
         ## fall back to the numeric gradient
         g <- numericGradient(fn, theta, ...)
         ## Check whether the gradient has enough rows.  This is the case
         ## if and only if loglik has enough rows, hence the error message
         ## about loglik.
         if( !is.matrix( g ) || nrow( g ) == 1 ) {
            stop( "if the gradients (argument 'grad') are not provided by the user,",
               " the BHHH method requires that the log-likelihood function",
               " (argument 'fn') returns a numeric vector,",
               " where each element must be the log-likelihood value corresponding",
               " to an individual (independent) observation" )
         }
         if( nrow( g ) < length( theta ) ) {
            stop( "the vector returned by the log-likelihood function",
               " (argument 'fn') must have at least as many elements",
               " as the number of parameters,",
               " where each element must be the log-likelihood value corresponding",
               " to an individual (independent) observation:\n",
               " currently, there are (is) ", length( theta ), " parameter(s)",
               " but the log likelihood function return only ", nrow( g ),
               " element(s)" )
         }
      }
      assign("gradVal", g, inherits=TRUE)
      return( g )
   }
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }
   a <- maxNR(fn, grad=gradient, hess=hess, start=start, iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
