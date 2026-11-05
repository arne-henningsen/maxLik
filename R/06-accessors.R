## Internal accessors for fitted 'maxLik' and 'maxim' objects.
##
## These functions deliberately avoid S3 dispatch.  Methods for subclasses
## of 'maxLik' or 'maxim' must not change how maxLik's own methods access
## the components of the fitted object.
##
## All accessors read raw components through .field()/unclass(), i.e. they
## do not dispatch on '$' or on any public generic.  Where the same
## quantity is stored under different names in fitted and summary objects
## (notably 'summary.maxLik'), the aliases are handled centrally here.

## Read one raw component without dispatching on the object's class
.field <- function(x, name) {
   unclass( x )[[name]]
}

## Return the value of the first non-NULL component among 'names'
.firstField <- function(x, names) {
   for( n in names ) {
      value <- .field( x, n )
      if( !is.null( value ) ) {
         return( value )
      }
   }
   NULL
}

## Parameter estimates (for 'summary.maxLik' this is the coefficient table)
.estimate <- function(x) {
   .field( x, "estimate" )
}

## Parameter names of a fitted object
.parameterNames <- function(x) {
   .field( x, "names" )
}

## Logical vector: which parameters are free (not fixed)
.activePar <- function(x) {
   if( !is.null( fixed <- .field( x, "fixed" ) ) ) {
      return( !fixed )
   }
   if( !is.null( activePar <- .field( x, "activePar" ) ) ) {
      return( activePar )
   }
   if( inherits( x, "maxLik" ) || inherits( x, "maxim" ) ) {
      return( rep( TRUE, length( .estimate( x ) ) ) )
   }
   ## fallback for arbitrary objects handled by activePar.default()
   NULL
}

## Hessian of a fitted object
.hessian <- function(x) {
   .field( x, "hessian" )
}

## Value of the objective function at the optimum
.maximum <- function(x) {
   .firstField( x, c( "maximum", "loglik" ) )
}

## Gradient at the optimum
.gradient <- function(x) {
   .field( x, "gradient" )
}

## Type of maximization
.maximType <- function(x) {
   .firstField( x, c( "type", "maximType" ) )
}

## Number of iterations
.nIter <- function(x) {
   .field( x, "iterations" )
}

## Return code
.returnCode <- function(x) {
   .firstField( x, c( "code", "returnCode" ) )
}

## Return message
.returnMessage <- function(x) {
   .firstField( x, c( "message", "returnMessage" ) )
}

## Optimization constraints
.constraints <- function(x) {
   .field( x, "constraints" )
}

## Observation-level gradients
.gradientObs <- function(x) {
   .field( x, "gradientObs" )
}

## Objective function used for optimization
.objectiveFn <- function(x) {
   .field( x, "objectiveFn" )
}

## Stored parameters at each iteration/epoch
.storedParameters <- function(x) {
   .field( x, "parameterStore" )
}

## Stored objective values at each iteration/epoch
.storedValues <- function(x) {
   .field( x, "valueStore" )
}

## Information about the last step
.lastStep <- function(x) {
   .field( x, "last.step" )
}

## Information about the last (unsuccessful) step in a summary object
.unsuccStep <- function(x) {
   .field( x, "unsucc.step" )
}

## Number of parameters, either all or only the free ones
.nParam <- function(x, free = FALSE) {
   if( free ) {
      sum( .activePar( x ) )
   } else {
      length( .estimate( x ) )
   }
}

## Number of free parameters of a summary object.  This is stored as
## 'NActivePar' in 'summary.maxLik' objects, but computed if missing.
.nActivePar <- function(x) {
   nActivePar <- .field( x, "NActivePar" )
   if( is.null( nActivePar ) ) {
      nActivePar <- sum( .activePar( x ) )
   }
   nActivePar
}

## Variance-covariance matrix of a fitted object
.vcov <- function(x, eigentol = 1e-12) {
   ## if exists $varcovar, take it
   ## NB: no maxLik function ever sets this component, and it is not
   ## documented.  It has behaved as a hook for user-supplied covariance
   ## matrices since at least version 0.5-2 (2008).  Whether it should be
   ## documented as such or deprecated is still undecided.
   if( !is.null( varcovar <- .field( x, "varcovar" ) ) ) {
      return( varcovar )
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

## Log-likelihood of a fitted or summary 'maxLik' object
.logLik <- function(x) {
   ll <- .maximum( x )
   attr( ll, "df" ) <- sum( .activePar( x ) )
   ll
}

## Akaike information criterion of a fitted 'maxLik' object
.AIC <- function(x, k = 2) {
   -2 * .logLik( x ) + k * .nParam( x, free = TRUE )
}

## Number of observations of a fitted 'maxLik' object
.nObs <- function(x) {
   gradientObs <- .gradientObs( x )
   if( is.null( gradientObs ) ) {
      stop( "cannot return the number of observations:",
         " please re-run 'maxLik' and",
         " provide a gradient function using argument 'grad' or",
         " (if no gradient function is specified) a log-likelihood function",
         " using argument 'logLik'",
         " that return the gradients or log-likelihood values, respectively,",
         " at each observation" )
   } else if( is.matrix( gradientObs ) ) {
      return( nrow( gradientObs ) )
   } else {
      stop( "internal error: component 'gradientObs' is not a matrix.",
         " Please contact the developers." )
   }
}
