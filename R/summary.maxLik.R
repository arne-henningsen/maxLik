print.summary.maxLik <- function( x,
      digits = max( 3L, getOption("digits") - 3L ), ... ) {
   
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(.maximType(x), ", ", .nIter(x), " iterations\n", sep="")
   cat("Return code ", .returnCode(x), ": ", .returnMessage(x), "\n", sep="")
   estimate <- .estimate( x )
   if(!is.null(estimate)) {
      cat("Log-Likelihood:", .maximum( x ), "\n")
      cat(.nActivePar( x ), " free parameters\n")
      cat("Estimates:\n")
      printCoefmat( estimate, digits = digits )
   }
   constraints <- .constraints( x )
   if(!is.null(constraints)) {
      cat("\nWarning: constrained likelihood estimation.",
          "Inference is probably wrong\n")
      cat("Constrained optimization based on", constraints$type,
          "\n")
      if(!is.null(constraints$code))
         cat("Return code:", constraints$code, "\n")
                           # note: this is missing for 'constrOptim'
      if(!is.null(constraints$message))
         cat(constraints$message, "\n")
                           # note: this is missing for 'constrOptim'
      cat(constraints$outer.iterations,
          " outer iterations, barrier value",
          constraints$barrier.value, "\n")
   }
   cat("--------------------------------------------\n")
}

summary.maxLik <- function(object, eigentol=1e-12,... ) {
   ## object      object of class "maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximType  : type of optimisation
   ## iterations : number of iterations
   ## returnCode : code of convergence
   ## returnMessage : message, description of the code
   ## loglik     : log-likelihood value at optimum
   ## estimate   : estimated parameter values, standard errors, t values
   ##              and p values at optimum
   ## fixed      : logical vector, which parameters are fixed
   ## NActivePar : number of free parameters
   ## constraints: information about constrained optimization
   ##
   if(!inherits(object, "maxLik"))
       stop("'summary.maxLik' called on a non-'maxLik' object")
   ## Use internal accessors, not the public generics.  This guarantees that
   ## methods for subclasses of 'maxLik' cannot affect the summary.
   estimate <- .estimate( object )
   activePar <- .activePar( object )
   code <- .returnCode( object )
   if((code < 100) & !is.null(estimate)) {
                           # in case of infinity at initial values, the coefs are not provided
       se <- .stdEr( object, eigentol = eigentol )
       t <- estimate / se
       p <- 2*pnorm( -abs( t))
       t[!activePar] <- NA
       p[!activePar] <- NA
       results <- cbind("Estimate" = estimate,
                        "Std. error" = se,
                        "t value" = t, "Pr(> t)" = p )
   } else {
     results <- NULL
   }
   summary <- list(maximType=.maximType( object ),
                   iterations=.nIter( object ),
                   returnCode=code,
                   returnMessage=.returnMessage( object ),
                   loglik=.maximum( object ),
                   estimate=results,
                   fixed=!activePar,
                   NActivePar=sum(activePar),
                   constraints=.constraints( object ))
   class(summary) <- "summary.maxLik"
   summary
}
