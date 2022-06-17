
print.maxim <- function(x, ...) {
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   cat("Value:", maxValue(x), " (", sum( activePar( x ) ), " free parameter(s))\n", sep = "")
   cat("Estimate(s):", coef(x), "\n")
}

setMethod("show", "maxim",
          function(object) print(object))
