
print.maxim <- function(x, ...) {
   cat(.maximType(x), ", ", .nIter(x), " iterations\n", sep="")
   cat("Return code ", .returnCode(x), ": ", .returnMessage(x), "\n", sep="")
   cat("Value:", .maximum(x), " (", sum( .activePar( x ) ), " free parameter(s))\n", sep = "")
   cat("Estimate(s):", .estimate(x), "\n")
}

## 'show' deliberately calls the S3 generic 'print()' rather than
## 'print.maxim()' directly: a subclass of 'maxim' that provides its own
## 'print' method should stay in control of how it is displayed.  The
## no-dispatch rule that applies to the extractor generics does not apply
## here, see the comments in 'inst/tinytest/test-subclass.R'.
setMethod("show", "maxim",
          function(object) print(object))
