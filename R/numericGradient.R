numericGradient <- function(f, t0, eps=1e-6, fixed, ...) {
   ## numeric gradient of a vector-valued function
   ## f           function, return Nval x 1 vector of values
   ## t0          NPar x 1 vector of parameters
   ## fixed       calculate the gradient based on these parameters only
   ## return:
   ## NvalxNPar matrix, gradient
   ## gradient along parameters which are not active are NA
   NPar <- length(t0)
   NVal <- length(f0 <- f(t0, ...))
   grad <- matrix(NA, NVal, NPar)
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
      grad[,i] <- (f(t2, ...) - f(t1, ...))/eps
   }
   return(grad)
}

