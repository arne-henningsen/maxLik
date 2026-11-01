print.maxLik <- function( x, ... ) {
   cat("Maximum Likelihood estimation\n")
   cat(.maximType(x), ", ", .nIter(x), " iterations\n", sep="")
   cat("Return code ", .returnCode(x), ": ", .returnMessage(x), "\n", sep="")
   estimate <- .estimate( x )
   if( !is.null( estimate ) ) {
      cat("Log-Likelihood:", .maximum( x ) )
      cat( " (", sum( .activePar( x ) ), " free parameter(s))\n", sep = "" )
      cat("Estimate(s):", estimate, "\n" )
   }
}

setMethod("show", "maxLik",
          function(object) print(object))
