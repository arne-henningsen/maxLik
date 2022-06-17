print.maxLik <- function( x, ... ) {
   cat("Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(coef(x))) {
      cat("Log-Likelihood:", maxValue(x) )
      cat( " (", sum( activePar( x ) ), " free parameter(s))\n", sep = "" )
      cat("Estimate(s):", coef(x), "\n" )
   }
}

setMethod("show", "maxLik",
          function(object) print(object))
