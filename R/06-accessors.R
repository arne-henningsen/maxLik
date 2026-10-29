## Internal accessors for fitted 'maxLik' and 'maxim' objects.
##
## These functions deliberately avoid S3 dispatch.  Methods for subclasses
## of 'maxLik' or 'maxim' must not change how maxLik's own methods access
## the components of the fitted object.

## Parameter estimates of a fitted object (not of a summary object)
.estimate <- function(x) {
   x$estimate
}

## Logical vector: which parameters are free (not fixed)
.activePar <- function(x) {
   if( !is.null( x$fixed ) ) {
      return( !x$fixed )
   }
   if( !is.null( x$activePar ) ) {
      return( x$activePar )
   }
   if( inherits( x, "maxLik" ) || inherits( x, "maxim" ) ) {
      return( rep( TRUE, length( x$estimate ) ) )
   }
   ## fallback for arbitrary objects handled by activePar.default()
   rep( TRUE, length( coef( x ) ) )
}

## Hessian of a fitted object
.hessian <- function(x) {
   x$hessian
}

## Number of parameters, either all or only the free ones
.nParam <- function(x, free = FALSE) {
   if( free ) {
      sum( .activePar( x ) )
   } else {
      length( .estimate( x ) )
   }
}

## Variance-covariance matrix of a fitted object
.vcov <- function(x, eigentol = 1e-12) {
   ## if exists $varcovar, take it
   if( !is.null( x$varcovar ) ) {
      return( x$varcovar )
   }

   ## otherwise invert the Hessian
   active <- .activePar( x )
   hess <- .hessian( x )
   if( !is.null( hess ) ) {
      ## if hessian exists, invert it
      hessActive <- hess[active, active, drop = FALSE]
      hessev <- abs( eigen( hessActive, symmetric = TRUE,
                           only.values = TRUE )$values )
      K <- .nParam( x )
      varcovar <- matrix( 0, K, K )
                           # this makes the fixed parameters to 0
      rownames( varcovar ) <- colnames( varcovar ) <- names( .estimate( x ) )
      if( min( hessev ) > ( eigentol * max( hessev ) ) ) {
      ## If the Hessian is not singular, fill in the free parameter values
         varcovar[active, active] <- solve( -hessActive )
         ## guarantee that the returned variance covariance matrix is symmetric
         varcovar <- ( varcovar + t( varcovar ) ) / 2
      } else {
      ## If singular, the free parameter values will be Inf
         varcovar[active, active] <- Inf
      }
      return( varcovar )
   } else {
      ## If hessian does not exist, return NA
      K <- .nParam( x )
      matrix( NA_real_, K, K )
   }
}

## Standard errors of a fitted object
.stdEr <- function(x, eigentol = 1e-12) {
   if( !is.null( vc <- .vcov( x, eigentol = eigentol ) ) ) {
      s <- sqrt( diag( vc ) )
      names( s ) <- names( .estimate( x ) )
      return( s )
   }
   ## if vcov is not working, return NULL
   NULL
}
