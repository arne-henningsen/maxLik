numericGradient <- function(f, t0, eps=1e-6, fixed,
                            ...) {
   ## numeric gradient of a vector-valued function
   ## f           function, return Nval x 1 vector of values
   ## t0          NPar x 1 vector of parameters
   ## fixed       calculate the gradient based on these parameters only
   ## return:
   ## NvalxNPar matrix, gradient
   ## gradient along parameters which are not active are NA
   NPar <- length(t0)
   nVal <- length(f0 <- f(t0, ...))
   grad <- matrix(NA, nVal, NPar)
   row.names(grad) <- names(f0)
   colnames(grad) <- names(t0)
   if(missing(fixed))
       fixed <- rep(FALSE, NPar)
   for(i in 1:NPar) {
      if(fixed[i])
          next
      t2 <- t1 <- t0
      t1[i] <- t0[i] - eps/2
      t2[i] <- t0[i] + eps/2
      ft1 <- f(t1, ...)
      ft2 <- f(t2, ...)
      ## give meaningful error message if the functions give vectors
      ## of different length at t1, t2
      if(length(ft1) == nVal & length(ft2) == nVal) {
         grad[,i] <- (ft2 - ft1)/eps
      }
      else {
         cat("Problem in numeric gradient\n")
         max.print <- 10
         if(length(ft1) != nVal) {
            cat("Function value at\n")
            print(t1[seq(length=min(max.print,length(t1)))])
            if(max.print < length(t1))
                cat("...\n")
            cat(":\n")
            print(ft1[seq(length=min(max.print,length(ft1)))])
            if(max.print < length(ft1))
                cat("...\n")
            cat("(length ", length(ft1), ") does not conform with ",
                "the original vector length ", nVal, "\n", sep="")
         }
         if(length(ft2) != nVal) {
            cat("Function value at\n")
            print(t2[seq(length=min(max.print,length(t2)))])
            if(max.print < length(t2))
                cat("...\n")
            cat(":\n")
            print(ft2[seq(length=min(max.print,length(ft2)))])
            if(max.print < length(ft2))
                cat("...\n")
            cat("(length ", length(ft1), ") does not conform with ",
                "the original vector length ", nVal, sep="")
         }
         cat("component", i, "will be set to NA\n")
         grad[,i] <- NA
      }
   }
   return(grad)
}

