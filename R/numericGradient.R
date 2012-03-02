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
         warnMsg <- "Problem in numeric gradient\n"
         max.print <- 10
         if(length(ft1) != nVal) {
            warnMsg <- c(warnMsg,"Function value at\n")
            warnMsg <- c(warnMsg,
                             format(t1[seq(length=min(max.print,length(t1)))]),
                             "\n")
            if(max.print < length(t1))
                warnMsg <- c(warnMsg, "...\n")
            warnMsg <- c(warnMsg, " =\n")
            warnMsg <- c(warnMsg,
                         format(paste(ft1[seq(length=min(max.print,length(ft1)))],
                                      collapse=" ")
                                ), "\n")
            if(max.print < length(ft1))
                warnMsg <- c(warnMsg, "...\n")
            warnMsg <- c(warnMsg, "(length ", length(ft1), ") does not conform with ",
                "the length at original value ", nVal, "\n", sep="")
         }
         if(length(ft2) != nVal) {
            warnMsg <- c(warnMsg, "Function value at\n")
            warnMsg <- c(warnMsg,
                         paste(format(t2[seq(length=min(max.print,length(t2)))]),
                               collapse=" "),
                         "\n")
            if(max.print < length(t2))
                warnMsg <- c(warnMsg,"...\n")
            warnMsg <- c(warnMsg, " =\n")
            warnMsg <- c(warnMsg,
                             paste(format(ft2[seq(length=min(max.print,length(ft2)))]),
                                   collapse=" "),
                             "\n")
            if(max.print < length(ft2))
                warnMsg <- c(warnMsg, "...\n")
            warnMsg <- c(warnMsg,
                         "(length ", length(ft2), ") does not conform with ",
                         "the length at original value ", nVal, "\n", sep="")
         }
         warnMsg <- c(warnMsg,"component ", i, " will be set to NA")
         warning(warnMsg)
         grad[,i] <- NA
      }
   }
   return(grad)
}

