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

## 'show' deliberately calls the S3 generic 'print()' rather than
## 'print.maxLik()' directly: a subclass of 'maxLik' that provides its own
## 'print' method should stay in control of how it is displayed.  The
## no-dispatch rule that applies to the extractor generics does not apply
## here, see the comments in 'inst/tinytest/test-subclass.R'.
setMethod("show", "maxLik",
          function(object) print(object)
          )
